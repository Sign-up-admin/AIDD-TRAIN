"""
测试前端停止功能 - 模拟前端行为
"""
import sys
import time
import requests
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'FLASH_DOCK-main'))
from services.compass_client import CompassClient

def test_stop_from_frontend():
    """模拟前端停止操作"""
    print("="*70)
    print("模拟前端停止操作测试")
    print("="*70)
    
    # 创建客户端
    client = CompassClient()
    
    # 获取任务列表
    print("\n1. 获取任务列表...")
    try:
        tasks = client.list_training_tasks()
        print(f"   找到 {len(tasks)} 个任务")
        
        # 查找运行中的任务
        running_tasks = [t for t in tasks if t.get('status') in ['running', 'initializing']]
        if not running_tasks:
            print("\n   没有运行中的任务")
            return 1
        
        print(f"\n   运行中的任务: {len(running_tasks)} 个")
        for task in running_tasks:
            print(f"   - {task.get('task_id')[:8]}... ({task.get('status')})")
        
        # 选择第一个运行中的任务
        selected_task = running_tasks[0]
        task_id = selected_task['task_id']
        print(f"\n2. 选择任务: {task_id}")
        
        # 模拟前端的状态检查
        print("\n3. 重新获取最新状态（模拟前端状态检查）...")
        try:
            current_task = client.get_training_task(task_id)
            current_status = current_task["status"]
            print(f"   当前状态: {current_status}")
            
            if current_status != "running" and current_status != "initializing":
                print(f"   ⚠️ 任务状态为 '{current_status}'，无法停止")
                return 1
            
            # 发送停止请求
            print("\n4. 发送停止请求...")
            start_time = time.time()
            try:
                result = client.stop_training_task(task_id)
                elapsed = time.time() - start_time
                print(f"   ✓ 停止请求已发送（耗时: {elapsed:.2f}秒）")
                print(f"   响应: {result}")
                
                # 轮询检查任务状态（模拟前端轮询）
                print("\n5. 轮询检查任务状态（最多5秒）...")
                max_poll_time = 5.0
                poll_interval = 0.5
                poll_elapsed = 0.0
                final_status = None
                
                while poll_elapsed < max_poll_time:
                    try:
                        updated_task = client.get_training_task(task_id)
                        final_status = updated_task.get('status')
                        print(f"   [{poll_elapsed:.1f}s] 状态: {final_status}")
                        
                        if final_status in ['cancelled', 'completed', 'failed']:
                            print(f"   ✓ 任务状态已更新为: {final_status}")
                            break
                        elif final_status not in ['running', 'initializing']:
                            print(f"   ⚠️ 任务状态: {final_status} (非预期状态)")
                            break
                    except Exception as poll_e:
                        print(f"   ⚠️ 轮询时出错: {poll_e}")
                        break
                    
                    time.sleep(poll_interval)
                    poll_elapsed += poll_interval
                
                if final_status and final_status in ['running', 'initializing']:
                    print(f"\n   ⚠️ 警告: 5秒后任务状态仍为 {final_status}")
                    print("   提示: 停止请求已发送，但任务可能仍在停止中")
                    print("   建议: 继续轮询或刷新页面查看最新状态")
                    
                    # 继续轮询更长时间
                    print("\n6. 继续轮询（额外30秒）...")
                    additional_poll_time = 30.0
                    additional_elapsed = 0.0
                    
                    while additional_elapsed < additional_poll_time:
                        try:
                            updated_task = client.get_training_task(task_id)
                            final_status = updated_task.get('status')
                            
                            if final_status != 'running' and final_status != 'initializing':
                                print(f"   [{additional_elapsed:.1f}s] 状态已更新为: {final_status}")
                                break
                        except Exception as poll_e:
                            print(f"   ⚠️ 轮询时出错: {poll_e}")
                            break
                        
                        time.sleep(poll_interval)
                        additional_elapsed += poll_interval
                    
                    if final_status in ['running', 'initializing']:
                        print(f"\n   ❌ 问题: 35秒后任务状态仍为 {final_status}")
                        print("   这可能表示停止功能存在问题")
                        return 1
                    else:
                        print(f"\n   ✓ 最终状态: {final_status}")
                        return 0
                else:
                    print(f"\n   ✓ 最终状态: {final_status}")
                    return 0
                    
            except Exception as e:
                elapsed = time.time() - start_time
                print(f"   ❌ 停止请求失败（耗时: {elapsed:.2f}秒）")
                print(f"   错误: {type(e).__name__}: {e}")
                import traceback
                traceback.print_exc()
                return 1
                
        except Exception as status_check_e:
            print(f"   ❌ 状态检查失败: {status_check_e}")
            return 1
            
    except Exception as e:
        print(f"\n❌ 错误: {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(test_stop_from_frontend())

