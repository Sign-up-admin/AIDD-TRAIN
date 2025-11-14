# COMPASS: Navigating the Complexities of Protein-Ligand Binding

[![Python Version](https://img.shields.io/badge/Python-3.12-blue.svg)](https://www.python.org/)
[![PyTorch](https://img.shields.io/badge/PyTorch-2.8-orange.svg)](https://pytorch.org/)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)

**COMPASS** is a deep learning project dedicated to accurately predicting protein-ligand binding affinities. It leverages a state-of-the-art Graph Neural Network (GNN), ViSNet, to learn from the intricate 3D geometry of molecular complexes, aiming to accelerate the process of drug discovery.

This project is built not just on a powerful model, but on a core philosophy of extreme data robustness and a highly efficient, mode-driven development workflow. Its name, COMPASS, reflects its ability to navigate the often-problematic landscape of real-world structural biology data, ensuring reliable and reproducible results.

---

## Key Features

- **State-of-the-Art Model**: Implements the ViSNet architecture for high-precision, geometry-aware predictions.
- **Robust Data Pipeline**: The cornerstone of COMPASS. The data processing pipeline is meticulously designed to handle common and obscure issues found in PDB data.
- **Automated Hardware Optimization**: A built-in tool to automatically find the best-performing configuration for your specific hardware, eliminating memory errors.
- **Mode-Driven Workflow**: A four-stage development process (`smoke_test`, `prototyping`, `validation`, `production`) that allows for seamless switching between quick checks, rapid experimentation, and full-scale training.
- **High-Performance Training**: Utilizes Automatic Mixed Precision (AMP) for significant speed-ups.
- **Resilient & Manageable**: Features robust checkpointing, graceful exit handling, and automated log/checkpoint organization.
- **Important ViSNet Constraint**: The ViSNet model requires that the number of hidden channels be divisible by the number of attention heads. This is a key consideration when configuring the model.

---

## Step 1: Setup and Installation

1.  **Clone the repository:**
    ```sh
    git clone https://github.com/Sign-up-admin/AIDD-TRAIN.git
    cd AIDD-TRAIN
    ```

2.  **Create a virtual environment (recommended):**
    ```sh
    python -m venv venv
    source venv/bin/activate  # On Windows, use `venv\Scripts\activate`
    ```

3.  **Install dependencies:**
    *This project relies on PyTorch and PyTorch Geometric. Please follow their official installation instructions for your specific CUDA version first.https://pytorch.org/get-started/locally/*
    ```sh
    # Example for CUDA 12.8 win11
    pip3 install torch torchvision --index-url https://download.pytorch.org/whl/cu128
    pip install torch_geometric
    pip install pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv -f https://data.pyg.org/whl/torch-2.8.0+cu128.html
    pip install rdkit-pypi biopython tqdm scipy
    ```

4.  **Prepare Data:**
    - Download the PDBbind dataset (v.2020 or other).
    - Open `compass/config.py` and update the `dataset_path` and `index_file` variables to point to your dataset location.

---

## Step 2: Automated Hardware Optimization (Highly Recommended)

This project contains an intelligent hardware optimizer, `hardware_optimizer.py`, designed to find the perfect model configuration for different stages of the development lifecycle.

### The Core Philosophy (see `docs/optimizer_philosophy.md` for details)

This optimizer was built upon a clear, hierarchical development philosophy defined by its architect, in collaboration with the Gemini agent. It recognizes that the "best" configuration is not a single setting, but a set of trade-offs tailored to the specific goal of each development phase. The optimizer's intelligence lies in its ability to weigh these trade-offs, using real-time performance estimation.

The optimizer targets three core stages, each with a unique goal:

1.  **`prototyping` (Soft Target: ~20 min/cycle)**
    *   **Philosophy**: **Time is the ruler.** This stage is for rapid trial and error. The configuration must be fast enough to allow developers to quickly test ideas. The optimizer targets a ~20 minute cycle time (for a fixed 450-batch run) but allows for a **20-minute flexibility window**. 
    *   **Strategy**: It searches its dedicated **small model space** to find the configuration with the **highest throughput (max batch size)** that fits within this flexible time budget. This ensures the fastest possible iteration speed without prematurely discarding a slightly slower but much more powerful configuration.

2.  **`validation` (Soft Target: ~90 min/cycle)**
    *   **Philosophy**: **Balance is the key.** This stage acts as the crucial bridge between a promising prototype and a full-scale production run. It must be close enough to production quality to give meaningful results, but fast enough to not halt the development flow. It serves to seriously validate the findings from the `prototyping` stage.
    *   **Strategy**: It targets a ~90 minute cycle time, also with a **20-minute flexibility window**. It searches its dedicated **large model space** for the configuration with the **highest throughput**, striking the perfect balance between speed and quality.

3.  **`production` (Goal: Time-unlimited)**
    *   **Philosophy**: **Quality is the ultimate goal.** Time is no longer the primary constraint. This stage is for building the best possible model that the hardware and data can support, ready for deployment.
    *   **Strategy**: It employs a **two-stage optimization**: first, it finds the highest-quality (largest) model that respects both data and hardware limits. Second, it squeezes all remaining performance out of the hardware by finding the maximum possible batch size for that single best model.

This structured approach ensures that from the earliest idea to the final deployment, there is a perfectly optimized configuration to support the task at hand.

### Running the Optimizer

1.  Open your terminal.
2.  Navigate to the project's root directory (`AIDD-TRAIN`).
3.  Run the optimizer module with the following command:

    ```sh
    python -m compass.optimizer
    ```

This command will optimize for all four modes (`production`, `validation`, `prototyping`, and `smoke_test`) in the correct order. The process may take some time.

The script will create or update a `hardware_profile.json` file in the project root. The main training script will automatically load the appropriate settings from this file based on the `DEVELOPMENT_MODE` you select in `compass/config.py`.

### (Optional) Optimizing for Specific Modes

If you wish to re-run the optimization for only specific modes, you can use the `--modes` argument:

```sh
# Example: Optimize only for production and validation
python -m compass.optimizer --modes production validation
```

---

## Step 3: Configuration and Daily Usage

Once the one-time setup and optimization are complete, your daily workflow is very simple.

1.  **Select Your Mode**: Open `compass/config.py` and set the `DEVELOPMENT_MODE` variable to one of the four modes: `'smoke_test'`, `'prototyping'`, `'validation'`, or `'production'`.

2.  **Run Training**: Execute the main script from your terminal.
    ```sh
    python -m compass
    ```

The script will automatically use the best settings for your chosen mode—either the optimized parameters from hardware_profile.json or the default settings if no optimization was run for that mode.

All logs and model checkpoints will be saved into a uniquely named directory (e.g., `checkpoints/visnet_prototyping_.../`).

---

## Understanding the Four-Stage Workflow

COMPASS implements a phased workflow to balance speed and rigor. You can switch between modes by changing a single variable in `compass/config.py`.

1.  **`smoke_test`**: *"Does the code run?"* A minimal check that runs in minutes.
2.  **`prototyping`**: *"Is my idea promising?"* A lightweight configuration for rapid experimentation.
3.  **`validation`**: *"How does my idea perform under realistic conditions?"* A medium-sized configuration for pre-production validation.
4.  **`production`**: *"What are the final, best-effort results?"* The full-scale configuration for generating final results.

---

## The COMPASS Philosophy: A Case Study in Data Robustness

This project was forged through a deep-dive debugging session to solve the sudden appearance of `NaN` (Not a Number) values during training. Instead of simply skipping problematic data, we developed a strategy of **"Pause and Autopsy"**.

This journey underscores the COMPASS philosophy: true progress in scientific machine learning comes not just from powerful architectures, but from a relentless commitment to understanding and purifying the data that fuels them.

---

## Security Configuration

### Production Deployment Security

For production deployments, the following security measures should be configured:

#### 1. CORS Configuration

Set the `CORS_ORIGINS` environment variable to specify allowed origins:

```bash
# Production example
export CORS_ORIGINS="https://yourdomain.com,https://api.yourdomain.com"

# Development (default)
# Uses: http://localhost:8501,http://127.0.0.1:8501,http://localhost:3000,http://127.0.0.1:3000
```

**Important**: Never use wildcard (`*`) origins in production. The system will reject wildcard origins for security.

#### 2. Authentication

Enable API key authentication for production:

```bash
# Enable authentication
export AUTH_ENABLED="true"

# Set API key (single key)
export API_KEY="your-secure-api-key-here"

# Or set multiple API keys (comma-separated) for key rotation
export API_KEYS="key1,key2,key3"

# Force authentication for critical endpoints in production
export FORCE_AUTH_CRITICAL="true"
```

**Critical endpoints** that require authentication in production:
- `/api/v1/training/tasks` - Training task management
- `/api/v1/data/upload` - Dataset uploads
- `/api/v1/data/datasets` - Dataset management
- `/api/v1/inference` - Model inference
- `/api/v1/models` - Model management

#### 3. Rate Limiting

Configure rate limits to prevent abuse:

```bash
# Default rate limit (requests per minute)
export RATE_LIMIT_DEFAULT="100"

# Training endpoints (more restrictive)
export RATE_LIMIT_TRAINING="10"
export RATE_LIMIT_TRAINING_WINDOW="60"

# Upload endpoints (very restrictive)
export RATE_LIMIT_UPLOAD="3"
export RATE_LIMIT_UPLOAD_WINDOW="60"

# Inference endpoints
export RATE_LIMIT_INFERENCE="20"
export RATE_LIMIT_INFERENCE_WINDOW="60"
```

#### 4. Environment Configuration

Set the environment to production:

```bash
export ENVIRONMENT="production"
```

This enables:
- Stricter security checks
- Enhanced error message sanitization
- Production-optimized logging
- Required authentication for critical endpoints

#### 5. Database Security

Configure database connection timeouts:

```bash
# Database connection timeout (seconds)
export DB_CONNECTION_TIMEOUT="10.0"

# Database busy timeout (milliseconds)
export DB_BUSY_TIMEOUT="5000"

# Database cache size (negative = KB)
export DB_CACHE_SIZE="-2000"
```

### Security Headers

The service automatically adds security headers to all responses:
- Content-Security-Policy (CSP) - Strict policy for API endpoints
- X-Frame-Options: DENY
- X-Content-Type-Options: nosniff
- X-XSS-Protection: 1; mode=block
- Referrer-Policy: strict-origin-when-cross-origin

### Input Validation

All user inputs are automatically sanitized to prevent:
- XSS (Cross-Site Scripting) attacks
- SQL injection (via parameterized queries)
- Path traversal attacks
- Command injection

### File Upload Security

File uploads are protected by:
- File type validation (only `.zip`, `.tar`, `.tar.gz` allowed)
- File size limits (configurable via `COMPASS_UPLOAD_MAX_SIZE`)
- Zip bomb detection
- Upload queue management (prevents resource exhaustion)

---

## Deployment Guide

### Production Deployment

#### 1. Prerequisites

- Python 3.12+
- CUDA-capable GPU (recommended)
- Sufficient disk space for datasets and checkpoints
- Network access for service registry (if using)

#### 2. Environment Setup

```bash
# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
pip install -r requirements_service.txt

# Install PyTorch with CUDA support
# See: https://pytorch.org/get-started/locally/
```

#### 3. Configuration

Create a `.env` file or set environment variables:

```bash
# Service configuration
export ENVIRONMENT="production"
export COMPASS_HOST="0.0.0.0"
export COMPASS_PORT="8080"

# Security
export AUTH_ENABLED="true"
export API_KEY="your-secure-api-key"
export CORS_ORIGINS="https://yourdomain.com"

# Resource limits
export COMPASS_MAX_WORKERS="4"
export MAX_CONCURRENT_UPLOADS="2"
export MAX_CONCURRENT_TASKS="4"

# Database
export REGISTRY_DB_PATH="./registry.db"
export DB_CONNECTION_TIMEOUT="10.0"

# Logging
export LOG_LEVEL="INFO"
export COMPASS_LOG_DIR="./logs"
```

#### 4. Start Services

```bash
# Start COMPASS service
python -m compass.service

# Or use the service startup script
python compass/service_main.py
```

#### 5. Health Checks

Monitor service health:

```bash
# Basic health check
curl http://localhost:8080/health

# Readiness check
curl http://localhost:8080/health/ready

# Metrics endpoint
curl http://localhost:8080/metrics
```

#### 6. Using with Reverse Proxy (Nginx)

Example Nginx configuration:

```nginx
server {
    listen 80;
    server_name yourdomain.com;

    location / {
        proxy_pass http://localhost:8080;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
    }
}
```

#### 7. Process Management

For production, use a process manager like `systemd` or `supervisord`:

**systemd example** (`/etc/systemd/system/compass.service`):

```ini
[Unit]
Description=COMPASS Service
After=network.target

[Service]
Type=simple
User=compass
WorkingDirectory=/path/to/AIDD-TRAIN
Environment="PATH=/path/to/venv/bin"
ExecStart=/path/to/venv/bin/python -m compass.service
Restart=always
RestartSec=10

[Install]
WantedBy=multi-user.target
```

### Troubleshooting

#### Common Issues

1. **Service won't start**
   - Check port availability: `netstat -an | grep 8080`
   - Verify environment variables are set correctly
   - Check logs in `logs/` directory

2. **Authentication failures**
   - Verify `AUTH_ENABLED` and `API_KEY` are set
   - Check API key format in request headers: `X-API-Key: your-key` or `Authorization: Bearer your-key`

3. **Rate limiting issues**
   - Check rate limit statistics: `curl http://localhost:8080/metrics`
   - Adjust rate limits via environment variables if needed

4. **Database connection errors**
   - Verify database file permissions
   - Check `DB_CONNECTION_TIMEOUT` setting
   - Ensure sufficient disk space

5. **File upload failures**
   - Check file size limits
   - Verify file type is allowed
   - Check upload queue capacity

#### Monitoring

Monitor service health and performance:

```bash
# Get metrics
curl http://localhost:8080/metrics

# Check rate limiting stats
curl http://localhost:8080/metrics | jq '.rate_limiting'

# View authentication failures (check logs)
tail -f logs/compass-service.log | grep "Authentication failed"
```

---

## Code Quality

This project maintains high code quality standards using automated tools:

### Quick Start

**Python Code Quality:**
```bash
# Install development dependencies
pip install -r requirements-dev.txt

# Run all code quality checks
python scripts/run_all_checks.py

# Auto-fix formatting issues
python scripts/run_all_checks.py --format

# Run tests
python scripts/run_tests.bat
```

**Frontend Code Quality:**

**Option 1: Using Docker (Recommended, no Node.js installation needed)**
```bash
# Install Docker Desktop first
# Then run checks using Docker
python scripts/check_frontend_docker.py
```

**Option 2: Using Local Node.js**
```bash
# Install Node.js dependencies (requires Node.js >= 14.0.0)
npm install

# Extract frontend code from Python files
python scripts/extract_frontend_code.py FLASH_DOCK-main

# Run frontend code checks
python scripts/check_frontend.py

# Or use npm scripts directly
npm run lint:all
npm run format
```

### Tools Used

**Python Code Quality:**
- **Black**: Code formatting (PEP 8)
- **Flake8**: Code style and complexity checking
- **Pylint**: Code quality analysis
- **MyPy**: Static type checking
- **Bandit**: Security vulnerability scanning
- **Pytest**: Unit testing and coverage

**Frontend Code Quality:**
- **ESLint**: JavaScript code linting
- **Stylelint**: CSS/SCSS code linting
- **Prettier**: Code formatting (HTML/CSS/JS)
- **HTMLHint**: HTML code quality checking

### Reports

All quality check reports are saved to `lint_reports/` directory. For detailed information, see:
- [Python Code Quality Guide](docs/CODE_QUALITY.md)
- [Frontend Code Quality Guide](docs/FRONTEND_CODE_QUALITY.md)

---

## License

This project is licensed under the GNU AGPLv3 License. See the `LICENSE` file for details.

## Acknowledgments

This project utilizes the PDBbind dataset. We gratefully acknowledge the creators and maintainers of this valuable resource.
