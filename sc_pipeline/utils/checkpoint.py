"""
断点管理模块

提供流程断点续传功能
"""

import os
import json
from typing import Optional, Dict
from datetime import datetime
import pandas as pd


class CheckpointManager:
    """
    断点管理器 - 支持从中断处继续运行

    功能：
    - 记录每个步骤的完成状态
    - 保存中间结果文件路径
    - 支持从任意断点恢复
    - 显示当前进度
    """

    def __init__(self, output_dir: str):
        """
        初始化断点管理器

        参数：
        ----------
        output_dir : str
            输出目录
        """
        self.output_dir = output_dir
        self.checkpoint_file = os.path.join(output_dir, '.checkpoint.json')
        self.checkpoints = self.load_checkpoints()

    def load_checkpoints(self) -> dict:
        """加载已有的检查点"""
        if os.path.exists(self.checkpoint_file):
            try:
                with open(self.checkpoint_file, 'r') as f:
                    return json.load(f)
            except:
                return {}
        return {}

    def save_checkpoint(self, step: str, **kwargs):
        """
        保存检查点

        参数：
        ----------
        step : str
            步骤名称
        **kwargs
            其他要保存的信息
        """
        self.checkpoints[step] = {
            'completed': True,
            'timestamp': datetime.now().isoformat(),
            **kwargs
        }
        with open(self.checkpoint_file, 'w') as f:
            json.dump(self.checkpoints, f, indent=2)
        print(f"   💾 断点已保存: {step}")

    def is_completed(self, step: str) -> bool:
        """
        检查某个步骤是否已完成

        参数：
        ----------
        step : str
            步骤名称

        返回：
        ----------
        bool
            是否已完成
        """
        return self.checkpoints.get(step, {}).get('completed', False)

    def get_checkpoint_data(self, step: str) -> dict:
        """
        获取检查点数据

        参数：
        ----------
        step : str
            步骤名称

        返回：
        ----------
        dict
            检查点数据
        """
        return self.checkpoints.get(step, {})

    def reset(self):
        """重置所有检查点（从头开始）"""
        self.checkpoints = {}
        if os.path.exists(self.checkpoint_file):
            os.remove(self.checkpoint_file)
        print("✅ 已重置所有检查点")

    def get_last_checkpoint(self) -> Optional[str]:
        """
        获取最后一个完成的步骤

        返回：
        ----------
        str or None
            最后完成的步骤名称
        """
        if not self.checkpoints:
            return None
        steps_order = ['step1_read', 'step2_qc', 'step3_integrate', 'step4_annotate']
        for step in reversed(steps_order):
            if self.is_completed(step):
                return step
        return None

    def print_status(self):
        """打印当前进度状态"""
        steps = [
            ('step1_read', '步骤1: 数据读取'),
            ('step2_qc', '步骤2: 质量控制'),
            ('step3_integrate', '步骤3: 数据整合'),
            ('step4_annotate', '步骤4: 细胞注释')
        ]

        print("\n" + "=" * 80)
        print("当前分析进度")
        print("=" * 80)
        for step_id, step_name in steps:
            if self.is_completed(step_id):
                data = self.get_checkpoint_data(step_id)
                timestamp = data.get('timestamp', '')
                file_path = data.get('file', '')
                print(f"✅ {step_name}")
                if timestamp:
                    print(f"   完成时间: {timestamp[:19]}")
                if file_path and os.path.exists(file_path):
                    size_mb = os.path.getsize(file_path) / 1024**2
                    print(f"   文件: {os.path.basename(file_path)} ({size_mb:.1f} MB)")
            else:
                print(f"⏳ {step_name}")
        print("=" * 80 + "\n")

    def get_progress_percentage(self) -> float:
        """
        获取进度百分比

        返回：
        ----------
        float
            进度百分比 (0-100)
        """
        total_steps = 4
        completed_steps = sum(1 for step in ['step1_read', 'step2_qc', 'step3_integrate', 'step4_annotate']
                            if self.is_completed(step))
        return (completed_steps / total_steps) * 100
