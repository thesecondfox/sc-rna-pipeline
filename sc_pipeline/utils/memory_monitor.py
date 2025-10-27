"""
内存监控模块

提供实时内存监控和内存使用跟踪功能
"""

import os
import psutil
import time
import threading
import logging
from typing import Optional, Callable
from datetime import datetime


class MemoryMonitor:
    """
    内存监控器，用于跟踪进程内存使用情况

    功能：
    - 实时监控内存使用
    - 设置内存阈值告警
    - 记录内存使用历史
    - 自动释放内存（可选）
    """

    def __init__(
        self,
        threshold_percent: float = 80.0,
        check_interval: float = 1.0,
        enable_logging: bool = True,
        log_file: Optional[str] = None
    ):
        """
        初始化内存监控器

        参数：
        ----------
        threshold_percent : float
            内存使用阈值百分比（默认 80%）
        check_interval : float
            检查间隔（秒，默认 1.0）
        enable_logging : bool
            是否启用日志记录（默认 True）
        log_file : str, optional
            日志文件路径（如果为 None，则输出到控制台）
        """
        self.threshold_percent = threshold_percent
        self.check_interval = check_interval
        self.enable_logging = enable_logging

        # 获取当前进程
        self.process = psutil.Process(os.getpid())

        # 内存使用历史
        self.memory_history = []

        # 监控状态
        self.is_monitoring = False
        self.monitor_thread = None

        # 回调函数
        self.threshold_callback: Optional[Callable] = None

        # 配置日志
        if enable_logging:
            self._setup_logging(log_file)

    def _setup_logging(self, log_file: Optional[str] = None):
        """配置日志系统"""
        self.logger = logging.getLogger('MemoryMonitor')
        self.logger.setLevel(logging.INFO)

        # 避免重复添加handler
        if not self.logger.handlers:
            if log_file:
                handler = logging.FileHandler(log_file)
            else:
                handler = logging.StreamHandler()

            formatter = logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
            )
            handler.setFormatter(formatter)
            self.logger.addHandler(handler)

    def get_current_memory(self) -> dict:
        """
        获取当前内存使用情况

        返回：
        ----------
        dict
            包含内存使用信息的字典
        """
        # 进程内存信息
        mem_info = self.process.memory_info()

        # 系统内存信息
        virtual_mem = psutil.virtual_memory()

        return {
            'timestamp': datetime.now(),
            'process_rss_mb': mem_info.rss / 1024 / 1024,  # 常驻内存 (MB)
            'process_vms_mb': mem_info.vms / 1024 / 1024,  # 虚拟内存 (MB)
            'process_percent': self.process.memory_percent(),  # 进程内存占比
            'system_total_mb': virtual_mem.total / 1024 / 1024,  # 系统总内存
            'system_available_mb': virtual_mem.available / 1024 / 1024,  # 系统可用内存
            'system_percent': virtual_mem.percent  # 系统内存使用率
        }

    def check_threshold(self) -> bool:
        """
        检查是否超过内存阈值

        返回：
        ----------
        bool
            如果超过阈值返回 True
        """
        mem_info = self.get_current_memory()

        if mem_info['system_percent'] >= self.threshold_percent:
            if self.enable_logging:
                self.logger.warning(
                    f"内存使用率达到 {mem_info['system_percent']:.1f}% "
                    f"(阈值: {self.threshold_percent}%)"
                )

            # 调用回调函数
            if self.threshold_callback:
                self.threshold_callback(mem_info)

            return True

        return False

    def _monitor_loop(self):
        """监控循环（在后台线程中运行）"""
        while self.is_monitoring:
            mem_info = self.get_current_memory()
            self.memory_history.append(mem_info)

            # 检查阈值
            self.check_threshold()

            time.sleep(self.check_interval)

    def start_monitoring(self):
        """启动后台监控"""
        if self.is_monitoring:
            self.logger.warning("监控已在运行中")
            return

        self.is_monitoring = True
        self.monitor_thread = threading.Thread(target=self._monitor_loop, daemon=True)
        self.monitor_thread.start()

        if self.enable_logging:
            self.logger.info("内存监控已启动")

    def stop_monitoring(self):
        """停止后台监控"""
        self.is_monitoring = False

        if self.monitor_thread:
            self.monitor_thread.join(timeout=self.check_interval * 2)

        if self.enable_logging:
            self.logger.info("内存监控已停止")

    def set_threshold_callback(self, callback: Callable):
        """
        设置阈值触发时的回调函数

        参数：
        ----------
        callback : Callable
            回调函数，接收内存信息字典作为参数
        """
        self.threshold_callback = callback

    def get_memory_summary(self) -> dict:
        """
        获取内存使用摘要

        返回：
        ----------
        dict
            内存使用统计摘要
        """
        if not self.memory_history:
            return self.get_current_memory()

        process_rss_values = [h['process_rss_mb'] for h in self.memory_history]
        system_percent_values = [h['system_percent'] for h in self.memory_history]

        import numpy as np

        return {
            'samples': len(self.memory_history),
            'current': self.get_current_memory(),
            'process_rss': {
                'min_mb': np.min(process_rss_values),
                'max_mb': np.max(process_rss_values),
                'mean_mb': np.mean(process_rss_values),
                'std_mb': np.std(process_rss_values)
            },
            'system_usage': {
                'min_percent': np.min(system_percent_values),
                'max_percent': np.max(system_percent_values),
                'mean_percent': np.mean(system_percent_values),
                'std_percent': np.std(system_percent_values)
            }
        }

    def print_summary(self):
        """打印内存使用摘要"""
        summary = self.get_memory_summary()

        print("\n" + "=" * 60)
        print("内存使用摘要")
        print("=" * 60)

        current = summary['current']
        print(f"\n当前状态:")
        print(f"  进程内存 (RSS): {current['process_rss_mb']:.1f} MB")
        print(f"  进程内存占比: {current['process_percent']:.1f}%")
        print(f"  系统内存使用率: {current['system_percent']:.1f}%")
        print(f"  系统可用内存: {current['system_available_mb']:.1f} MB")

        if summary['samples'] > 1:
            print(f"\n历史统计 (共 {summary['samples']} 个样本):")
            print(f"  进程内存 (RSS):")
            print(f"    最小值: {summary['process_rss']['min_mb']:.1f} MB")
            print(f"    最大值: {summary['process_rss']['max_mb']:.1f} MB")
            print(f"    平均值: {summary['process_rss']['mean_mb']:.1f} MB")
            print(f"  系统内存使用率:")
            print(f"    最小值: {summary['system_usage']['min_percent']:.1f}%")
            print(f"    最大值: {summary['system_usage']['max_percent']:.1f}%")
            print(f"    平均值: {summary['system_usage']['mean_percent']:.1f}%")

        print("=" * 60 + "\n")

    def clear_history(self):
        """清除内存使用历史"""
        self.memory_history = []

        if self.enable_logging:
            self.logger.info("内存历史已清除")


def log_memory(prefix: str = ""):
    """
    快速记录当前内存使用（装饰器或独立函数）

    参数：
    ----------
    prefix : str
        日志前缀
    """
    monitor = MemoryMonitor(enable_logging=False)
    mem_info = monitor.get_current_memory()

    print(f"{prefix}内存使用: "
          f"进程={mem_info['process_rss_mb']:.1f}MB "
          f"({mem_info['process_percent']:.1f}%), "
          f"系统={mem_info['system_percent']:.1f}%")


def memory_profiler(func):
    """
    内存分析装饰器

    用法：
    @memory_profiler
    def my_function():
        ...
    """
    def wrapper(*args, **kwargs):
        monitor = MemoryMonitor(enable_logging=False)

        # 记录开始内存
        start_mem = monitor.get_current_memory()
        print(f"\n[{func.__name__}] 开始执行")
        print(f"  初始内存: {start_mem['process_rss_mb']:.1f} MB")

        # 执行函数
        start_time = time.time()
        result = func(*args, **kwargs)
        elapsed_time = time.time() - start_time

        # 记录结束内存
        end_mem = monitor.get_current_memory()
        mem_increase = end_mem['process_rss_mb'] - start_mem['process_rss_mb']

        print(f"[{func.__name__}] 执行完成")
        print(f"  最终内存: {end_mem['process_rss_mb']:.1f} MB")
        print(f"  内存增量: {mem_increase:+.1f} MB")
        print(f"  执行时间: {elapsed_time:.2f} 秒\n")

        return result

    return wrapper
