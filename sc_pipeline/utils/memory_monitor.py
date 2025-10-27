"""
内存监控模块

提供实时内存监控和内存使用跟踪功能
"""

import os
import psutil
import pandas as pd
import matplotlib.pyplot as plt
from typing import Optional, Callable, List
from datetime import datetime


class MemoryMonitor:
    """
    内存监控器，用于跟踪进程内存使用情况

    功能：
    - 实时监控内存使用
    - 记录内存使用历史
    - 生成内存使用报告和图表
    """

    def __init__(self):
        """初始化内存监控器"""
        self.process = psutil.Process()
        self.checkpoints = []

    def get_memory_usage(self) -> float:
        """
        获取当前内存使用（GB）

        返回：
        ----------
        float
            内存使用量（GB）
        """
        return self.process.memory_info().rss / 1024**3

    def checkpoint(self, step_name: str):
        """
        记录检查点

        参数：
        ----------
        step_name : str
            步骤名称
        """
        mem_gb = self.get_memory_usage()
        self.checkpoints.append({'step': step_name, 'memory_gb': mem_gb})
        print(f"   📊 内存使用: {mem_gb:.2f} GB ({step_name})")

    def get_summary(self) -> pd.DataFrame:
        """
        获取内存使用摘要

        返回：
        ----------
        pd.DataFrame
            内存使用历史数据框
        """
        df = pd.DataFrame(self.checkpoints)
        if len(df) > 0:
            df['memory_increase_gb'] = df['memory_gb'].diff()
        return df

    def plot_memory_usage(self, output_path: str):
        """
        绘制内存使用趋势图

        参数：
        ----------
        output_path : str
            输出文件路径
        """
        df = self.get_summary()
        if len(df) == 0:
            return

        fig, ax = plt.subplots(figsize=(12, 6))
        ax.plot(range(len(df)), df['memory_gb'], marker='o', linewidth=2, markersize=8)
        ax.set_xticks(range(len(df)))
        ax.set_xticklabels(df['step'], rotation=45, ha='right')
        ax.set_ylabel('Memory Usage (GB)', fontsize=12)
        ax.set_xlabel('Pipeline Step', fontsize=12)
        ax.set_title('Memory Usage Throughout Pipeline', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"   💾 内存使用图已保存: {output_path}")


def log_memory(prefix: str = ""):
    """
    快速记录当前内存使用

    参数：
    ----------
    prefix : str
        日志前缀
    """
    monitor = MemoryMonitor()
    mem_gb = monitor.get_memory_usage()
    print(f"{prefix}内存使用: {mem_gb:.2f} GB")


def memory_profiler(func):
    """
    内存分析装饰器

    用法：
    @memory_profiler
    def my_function():
        ...
    """
    def wrapper(*args, **kwargs):
        import time
        monitor = MemoryMonitor()

        # 记录开始内存
        start_mem = monitor.get_memory_usage()
        print(f"\n[{func.__name__}] 开始执行")
        print(f"  初始内存: {start_mem:.2f} GB")

        # 执行函数
        start_time = time.time()
        result = func(*args, **kwargs)
        elapsed_time = time.time() - start_time

        # 记录结束内存
        end_mem = monitor.get_memory_usage()
        mem_increase = end_mem - start_mem

        print(f"[{func.__name__}] 执行完成")
        print(f"  最终内存: {end_mem:.2f} GB")
        print(f"  内存增量: {mem_increase:+.2f} GB")
        print(f"  执行时间: {elapsed_time:.2f} 秒\n")

        return result

    return wrapper
