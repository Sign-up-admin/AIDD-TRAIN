"""检查现有的训练任务"""
import requests
import sys

try:
    response = requests.get("http://localhost:8080/api/v1/training/tasks", timeout=10)
    if response.status_code == 200:
        data = response.json()
        # Handle different response formats
        if isinstance(data, list):
            tasks = data
        elif isinstance(data, dict) and 'tasks' in data:
            tasks = data['tasks']
        else:
            tasks = [data] if data else []
        
        print(f"Found {len(tasks)} tasks:")
        print()
        
        running_tasks = []
        for task in tasks:
            if isinstance(task, str):
                # If task is just a string (task_id), skip
                continue
            task_id = task.get('task_id', 'unknown') if isinstance(task, dict) else str(task)
            status = task.get('status', 'unknown') if isinstance(task, dict) else 'unknown'
            description = task.get('description', '') if isinstance(task, dict) else ''
            print(f"  Task ID: {task_id[:8] if len(task_id) > 8 else task_id}...")
            print(f"    Status: {status}")
            if description:
                print(f"    Description: {description}")
            print()
            
            if status in ['running', 'initializing']:
                running_tasks.append((task_id, status))
        
        if running_tasks:
            print(f"\nFound {len(running_tasks)} running/initializing tasks:")
            for task_id, status in running_tasks:
                print(f"  - {task_id} ({status})")
            print("\nYou can test stop functionality with these tasks!")
            print("Run: python debug_training_stop.py")
        else:
            print("\nNo running tasks found.")
            print("Cannot test stop functionality without a running task.")
            print("\nNote: CPU usage is too high (100%) to create new tasks.")
            print("Please wait for CPU usage to decrease or close other applications.")
    else:
        print(f"Failed to get tasks: {response.status_code}")
        print(response.text)
except Exception as e:
    print(f"Error: {e}")
    sys.exit(1)

