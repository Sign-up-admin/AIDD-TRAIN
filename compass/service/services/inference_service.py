"""
Inference service implementation.
"""
import os
import time
import logging
import threading
from typing import Dict, Optional, List, OrderedDict
from pathlib import Path
from collections import OrderedDict as OrderedDictType

import torch
from torch_geometric.data import Data, Batch

from compass.service.models.model import InferenceRequest, InferenceResponse, BatchInferenceResponse
from compass.service.services.model_service import ModelService
from compass.training.model import ViSNetPDB
from compass.training.checkpoint import load_latest_checkpoint_file
from compass.data.processing import process_item
from compass.data.loader.paths import get_pdb_info, get_data_paths
from compass.config import get_config
from compass.service.config import SERVICE_CONFIG

logger = logging.getLogger(__name__)


class LRUModelCache:
    """LRU cache for loaded models."""
    
    def __init__(self, max_size: int = 3):
        """
        Initialize LRU cache.
        
        Args:
            max_size: Maximum number of models to cache
        """
        self.max_size = max_size
        self.cache: OrderedDictType[str, torch.nn.Module] = OrderedDictType()
        self.lock = threading.Lock()
        self.access_times: Dict[str, float] = {}
    
    def get(self, key: str) -> Optional[torch.nn.Module]:
        """Get model from cache, moving it to end (most recently used)."""
        with self.lock:
            if key in self.cache:
                # Move to end (mark as recently used)
                model = self.cache.pop(key)
                self.cache[key] = model
                self.access_times[key] = time.time()
                return model
            return None
    
    def put(self, key: str, model: torch.nn.Module):
        """Put model in cache, evicting least recently used if needed."""
        with self.lock:
            if key in self.cache:
                # Update existing
                self.cache.move_to_end(key)
            else:
                # Check if we need to evict
                if len(self.cache) >= self.max_size:
                    # Remove least recently used
                    lru_key = next(iter(self.cache))
                    lru_model = self.cache.pop(lru_key)
                    del self.access_times[lru_key]
                    
                    # Clear GPU memory if on CUDA
                    if torch.cuda.is_available():
                        del lru_model
                        torch.cuda.empty_cache()
                    
                    logger.info(f"Evicted model {lru_key} from cache (LRU)")
                
                # Add new model
                self.cache[key] = model
                self.access_times[key] = time.time()
                logger.debug(f"Cached model {key} (cache size: {len(self.cache)}/{self.max_size})")
    
    def remove(self, key: str) -> bool:
        """Remove model from cache."""
        with self.lock:
            if key in self.cache:
                model = self.cache.pop(key)
                del self.access_times[key]
                
                # Clear GPU memory
                if torch.cuda.is_available():
                    del model
                    torch.cuda.empty_cache()
                
                logger.info(f"Removed model {key} from cache")
                return True
            return False
    
    def clear(self):
        """Clear all cached models."""
        with self.lock:
            for model in self.cache.values():
                if torch.cuda.is_available():
                    del model
            self.cache.clear()
            self.access_times.clear()
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
            logger.info("Cleared all cached models")
    
    def size(self) -> int:
        """Get current cache size."""
        with self.lock:
            return len(self.cache)
    
    def get_cached_keys(self) -> List[str]:
        """Get list of cached model IDs."""
        with self.lock:
            return list(self.cache.keys())


class InferenceService:
    """Service for model inference."""
    
    def __init__(self, model_service: ModelService):
        """
        Initialize inference service.
        
        Args:
            model_service: Model service instance
        """
        self.model_service = model_service
        max_cache_size = int(os.getenv('MODEL_CACHE_SIZE', '3'))
        self.model_cache = LRUModelCache(max_size=max_cache_size)
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        logger.info(f"Inference service initialized with device: {self.device}, cache size: {max_cache_size}")
    
    def load_model(self, model_id: Optional[str] = None) -> torch.nn.Module:
        """
        Load a model for inference.
        
        This method uses LRU caching to avoid reloading models. If the model
        is already in cache, it returns the cached version. Otherwise, it loads
        the model from disk and adds it to the cache.
        
        Args:
            model_id: Optional model ID to load. If None, loads the latest model.
            
        Returns:
            torch.nn.Module: Loaded model ready for inference
            
        Raises:
            FileNotFoundError: If model file doesn't exist
            RuntimeError: If model loading fails
        """
        # Get model info to determine cache key
        if model_id:
            model_info = self.model_service.get_model(model_id)
            cache_key = model_id
        else:
            model_info = self.model_service.get_latest_model()
            if model_info:
                cache_key = model_info.model_id
            else:
                cache_key = None
        
        if not model_info:
            raise ValueError("No model available for inference")
        
        # Check cache first
        cached_model = self.model_cache.get(cache_key)
        if cached_model is not None:
            logger.debug(f"Using cached model {cache_key}")
            return cached_model
        
        # Load from disk
        logger.info(f"Loading model {cache_key} from disk")
        
        # Load checkpoint
        checkpoint_dir = Path(model_info.checkpoint_path).parent
        checkpoint = load_latest_checkpoint_file(checkpoint_dir, self.device)
        
        if not checkpoint:
            raise ValueError(f"Could not load checkpoint from {checkpoint_dir}")
        
        # Get model config from checkpoint or use defaults
        config = get_config()
        
        # Create model
        model = ViSNetPDB(
            hidden_channels=config.get('visnet_hidden_channels', 128),
            num_layers=config.get('visnet_num_layers', 6),
            num_rbf=config.get('visnet_num_rbf', 64),
            cutoff=config.get('visnet_cutoff', 8.0),
            max_num_neighbors=config.get('max_num_neighbors', 32),
            lmax=config.get('visnet_lmax', 1),
            vecnorm_type=config.get('visnet_vecnorm_type', 'max_min')
        ).to(self.device)
        
        # Load weights
        model.load_state_dict(checkpoint['model_state_dict'])
        model.eval()
        
        # Cache model (LRU will evict if needed)
        self.model_cache.put(cache_key, model)
        
        logger.info(f"Loaded and cached model {cache_key}")
        return model
    
    def unload_model(self, model_id: str) -> bool:
        """
        Unload a model from cache.
        
        Args:
            model_id: Model ID to unload
            
        Returns:
            bool: True if model was unloaded
        """
        return self.model_cache.remove(model_id)
    
    def clear_cache(self):
        """Clear all cached models."""
        self.model_cache.clear()
    
    def get_cache_info(self) -> Dict:
        """Get cache information."""
        return {
            "size": self.model_cache.size(),
            "max_size": self.model_cache.max_size,
            "cached_models": self.model_cache.get_cached_keys()
        }
    
    def predict(self, request: InferenceRequest) -> InferenceResponse:
        """
        Perform single prediction.
        
        Args:
            request: Inference request
            
        Returns:
            InferenceResponse: Prediction result
        """
        start_time = time.time()
        
        try:
            # Load model
            model = self.load_model(request.model_id)
            
            # Process data
            data = self._process_input(request.protein_path, request.ligand_path)
            if data is None:
                raise ValueError("Failed to process input data")
            
            # Run inference
            with torch.no_grad():
                data = data.to(self.device)
                output = model(data)
                binding_affinity = output.item()
            
            inference_time = time.time() - start_time
            
            model_id = request.model_id or self.model_service.get_latest_model().model_id
            
            return InferenceResponse(
                binding_affinity=binding_affinity,
                model_id=model_id,
                inference_time=inference_time
            )
        except Exception as e:
            logger.error(f"Inference failed: {e}", exc_info=True)
            raise
    
    def batch_predict(self, requests: List[InferenceRequest]) -> BatchInferenceResponse:
        """
        Perform batch predictions.
        
        Args:
            requests: List of inference requests
            
        Returns:
            BatchInferenceResponse: Batch prediction results
        """
        if not requests:
            raise ValueError("Empty request list")
        
        results = []
        success_count = 0
        failed_count = 0
        
        # Pre-load model before processing to fail fast if model unavailable
        model_id = requests[0].model_id if requests else None
        try:
            model = self.load_model(model_id)
            logger.info(f"Loaded model {model_id or 'latest'} for batch prediction")
        except Exception as e:
            logger.error(f"Failed to load model for batch prediction: {e}", exc_info=True)
            # Return error for all requests
            for request in requests:
                results.append({
                    'protein_path': request.protein_path,
                    'ligand_path': request.ligand_path,
                    'error': 'Failed to load model for batch prediction',
                    'success': False
                })
                failed_count += 1
            
            return BatchInferenceResponse(
                results=results,
                total_count=len(requests),
                success_count=0,
                failed_count=failed_count
            )
        
        # Process each request using pre-loaded model
        for idx, request in enumerate(requests):
            try:
                # Use pre-loaded model instead of calling predict() which would reload
                data = self._process_input(request.protein_path, request.ligand_path)
                if data is None:
                    raise ValueError("Failed to process input data")
                
                # Run inference with pre-loaded model
                with torch.no_grad():
                    data = data.to(self.device)
                    output = model(data)
                    binding_affinity = output.item()
                
                results.append({
                    'protein_path': request.protein_path,
                    'ligand_path': request.ligand_path,
                    'binding_affinity': binding_affinity,
                    'success': True
                })
                success_count += 1
            except Exception as e:
                logger.warning(f"Failed to process request {idx + 1}/{len(requests)}: {e}")
                results.append({
                    'protein_path': request.protein_path,
                    'ligand_path': request.ligand_path,
                    'error': 'Failed to process request',
                    'success': False
                })
                failed_count += 1
        
        logger.info(f"Batch prediction completed: {success_count} succeeded, {failed_count} failed out of {len(requests)}")
        
        return BatchInferenceResponse(
            results=results,
            total_count=len(requests),
            success_count=success_count,
            failed_count=failed_count
        )
    
    def _process_input(self, protein_path: str, ligand_path: str) -> Optional[Data]:
        """
        Process input files into graph data.
        
        Args:
            protein_path: Path to protein PDB file
            ligand_path: Path to ligand SDF file
            
        Returns:
            Optional[Data]: Processed data or None if failed
        """
        try:
            # Create item dict for processing
            item = {
                'pdb_code': Path(protein_path).stem,
                'protein_path': protein_path,
                'ligand_path': ligand_path,
                'binding_data': '0.0'  # Dummy value for inference
            }
            
            # Process item
            data = process_item(item)
            return data
        except Exception as e:
            logger.error(f"Failed to process input: {e}", exc_info=True)
            return None

