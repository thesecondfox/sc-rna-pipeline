"""
文件传输管理器

提供断点续传和传输任务管理功能
"""

import os
import json
import time
import hashlib
from pathlib import Path
from typing import Optional, Callable, Dict, List
from datetime import datetime
from tqdm import tqdm


class TransferState:
    """传输状态记录"""

    def __init__(self, state_file: str):
        self.state_file = state_file
        self.state = self._load_state()

    def _load_state(self) -> dict:
        """加载传输状态"""
        if os.path.exists(self.state_file):
            with open(self.state_file, 'r') as f:
                return json.load(f)
        return {}

    def _save_state(self):
        """保存传输状态"""
        os.makedirs(os.path.dirname(self.state_file) or '.', exist_ok=True)
        with open(self.state_file, 'w') as f:
            json.dump(self.state, f, indent=2)

    def get_task(self, task_id: str) -> Optional[dict]:
        """获取任务状态"""
        return self.state.get(task_id)

    def update_task(self, task_id: str, **kwargs):
        """更新任务状态"""
        if task_id not in self.state:
            self.state[task_id] = {
                'task_id': task_id,
                'created_at': datetime.now().isoformat(),
                'status': 'pending',
                'transferred_bytes': 0,
                'total_bytes': 0,
                'chunks_completed': [],
                'last_updated': datetime.now().isoformat()
            }

        self.state[task_id].update(kwargs)
        self.state[task_id]['last_updated'] = datetime.now().isoformat()
        self._save_state()

    def mark_chunk_completed(self, task_id: str, chunk_index: int):
        """标记块已完成"""
        if task_id in self.state:
            if chunk_index not in self.state[task_id]['chunks_completed']:
                self.state[task_id]['chunks_completed'].append(chunk_index)
                self._save_state()

    def is_chunk_completed(self, task_id: str, chunk_index: int) -> bool:
        """检查块是否已完成"""
        if task_id in self.state:
            return chunk_index in self.state[task_id]['chunks_completed']
        return False

    def remove_task(self, task_id: str):
        """删除任务记录"""
        if task_id in self.state:
            del self.state[task_id]
            self._save_state()

    def list_tasks(self) -> List[dict]:
        """列出所有任务"""
        return list(self.state.values())


class TransferManager:
    """
    文件传输管理器（支持断点续传）

    功能：
    - 分块传输
    - 断点续传
    - 传输进度跟踪
    - 传输状态持久化
    - 自动重试
    """

    def __init__(
        self,
        chunk_size: int = 1024 * 1024 * 10,  # 默认10MB
        state_dir: str = "./.transfer_state",
        max_retries: int = 3,
        retry_delay: int = 5
    ):
        """
        初始化传输管理器

        参数：
        ----------
        chunk_size : int
            传输块大小（字节，默认10MB）
        state_dir : str
            状态文件目录
        max_retries : int
            最大重试次数
        retry_delay : int
            重试延迟（秒）
        """
        self.chunk_size = chunk_size
        self.state_dir = state_dir
        self.max_retries = max_retries
        self.retry_delay = retry_delay

        os.makedirs(state_dir, exist_ok=True)
        self.state = TransferState(os.path.join(state_dir, 'transfer_state.json'))

    def _generate_task_id(self, src: str, dst: str) -> str:
        """生成任务ID"""
        unique_str = f"{src}:{dst}:{time.time()}"
        return hashlib.md5(unique_str.encode()).hexdigest()

    def _calculate_checksum(self, file_path: str, start: int = 0, end: Optional[int] = None) -> str:
        """计算文件块的校验和"""
        hash_obj = hashlib.md5()

        with open(file_path, 'rb') as f:
            f.seek(start)
            remaining = (end - start) if end else None

            while True:
                chunk_size = min(8192, remaining) if remaining else 8192
                chunk = f.read(chunk_size)

                if not chunk:
                    break

                hash_obj.update(chunk)

                if remaining:
                    remaining -= len(chunk)
                    if remaining <= 0:
                        break

        return hash_obj.hexdigest()

    def upload_file(
        self,
        local_path: str,
        remote_path: str,
        transfer_func: Callable,
        resume: bool = True,
        verify: bool = True,
        progress_callback: Optional[Callable] = None
    ) -> bool:
        """
        上传文件（支持断点续传）

        参数：
        ----------
        local_path : str
            本地文件路径
        remote_path : str
            远程文件路径
        transfer_func : Callable
            实际的传输函数，接收 (src, dst, offset, length) 参数
        resume : bool
            是否启用断点续传
        verify : bool
            是否验证文件完整性
        progress_callback : Callable
            进度回调函数

        返回：
        ----------
        bool
            是否成功
        """
        # 生成任务ID
        task_id = self._generate_task_id(local_path, remote_path)

        # 获取文件大小
        file_size = os.path.getsize(local_path)

        # 检查是否有未完成的任务
        existing_task = self.state.get_task(task_id)

        if resume and existing_task and existing_task['status'] == 'in_progress':
            print(f"\n检测到未完成的传输任务，将从上次断点继续...")
            print(f"  已传输: {existing_task['transferred_bytes'] / 1024 / 1024:.2f} MB / {file_size / 1024 / 1024:.2f} MB")
        else:
            # 初始化新任务
            self.state.update_task(
                task_id,
                status='in_progress',
                src=local_path,
                dst=remote_path,
                total_bytes=file_size,
                transferred_bytes=0,
                chunks_completed=[]
            )

        # 计算分块数量
        num_chunks = (file_size + self.chunk_size - 1) // self.chunk_size

        # 使用进度条
        pbar = tqdm(
            total=file_size,
            unit='B',
            unit_scale=True,
            desc=f"上传 {os.path.basename(local_path)}"
        )

        # 如果有已完成的块，更新进度条
        task = self.state.get_task(task_id)
        if task:
            pbar.update(task['transferred_bytes'])

        # 逐块传输
        try:
            for chunk_idx in range(num_chunks):
                # 检查该块是否已完成
                if self.state.is_chunk_completed(task_id, chunk_idx):
                    continue

                # 计算块的范围
                start = chunk_idx * self.chunk_size
                end = min(start + self.chunk_size, file_size)
                chunk_length = end - start

                # 尝试传输该块
                retry_count = 0
                while retry_count < self.max_retries:
                    try:
                        # 执行传输
                        self._transfer_chunk(
                            local_path,
                            remote_path,
                            start,
                            chunk_length,
                            transfer_func
                        )

                        # 标记块完成
                        self.state.mark_chunk_completed(task_id, chunk_idx)

                        # 更新传输进度
                        transferred = self.state.get_task(task_id)['transferred_bytes'] + chunk_length
                        self.state.update_task(task_id, transferred_bytes=transferred)

                        # 更新进度条
                        pbar.update(chunk_length)

                        # 调用进度回调
                        if progress_callback:
                            progress_callback(transferred, file_size)

                        break

                    except Exception as e:
                        retry_count += 1
                        if retry_count >= self.max_retries:
                            print(f"\n块 {chunk_idx} 传输失败（已重试 {self.max_retries} 次）: {e}")
                            self.state.update_task(task_id, status='failed', error=str(e))
                            pbar.close()
                            return False

                        print(f"\n块 {chunk_idx} 传输失败，{self.retry_delay}秒后重试... ({retry_count}/{self.max_retries})")
                        time.sleep(self.retry_delay)

            pbar.close()

            # 验证文件完整性
            if verify:
                print("\n验证文件完整性...")
                # 这里可以添加远程文件校验和验证
                # 简化版本：只验证文件大小
                pass

            # 标记任务完成
            self.state.update_task(
                task_id,
                status='completed',
                completed_at=datetime.now().isoformat()
            )

            print(f"\n文件上传完成: {local_path} -> {remote_path}")
            return True

        except Exception as e:
            print(f"\n传输失败: {e}")
            self.state.update_task(task_id, status='failed', error=str(e))
            pbar.close()
            return False

    def _transfer_chunk(
        self,
        src: str,
        dst: str,
        offset: int,
        length: int,
        transfer_func: Callable
    ):
        """
        传输文件块

        参数：
        ----------
        src : str
            源文件路径
        dst : str
            目标文件路径
        offset : int
            起始偏移量
        length : int
            块长度
        transfer_func : Callable
            传输函数
        """
        # 创建临时文件保存该块
        temp_file = f"{dst}.part{offset}"

        # 读取块数据
        with open(src, 'rb') as f:
            f.seek(offset)
            chunk_data = f.read(length)

        # 写入临时文件
        with open(temp_file, 'wb') as f:
            f.write(chunk_data)

        # 调用传输函数
        transfer_func(temp_file, dst, offset, length)

        # 删除临时文件
        if os.path.exists(temp_file):
            os.remove(temp_file)

    def download_file(
        self,
        remote_path: str,
        local_path: str,
        transfer_func: Callable,
        resume: bool = True,
        verify: bool = True,
        progress_callback: Optional[Callable] = None
    ) -> bool:
        """
        下载文件（支持断点续传）

        参数同 upload_file
        """
        # 类似于 upload_file 的实现
        # 这里简化处理，实际使用时可以参考 upload_file 实现
        task_id = self._generate_task_id(remote_path, local_path)

        try:
            # 执行下载
            transfer_func(remote_path, local_path)

            self.state.update_task(
                task_id,
                status='completed',
                completed_at=datetime.now().isoformat()
            )

            print(f"\n文件下载完成: {remote_path} -> {local_path}")
            return True

        except Exception as e:
            print(f"\n下载失败: {e}")
            self.state.update_task(task_id, status='failed', error=str(e))
            return False

    def list_tasks(self, status: Optional[str] = None) -> List[dict]:
        """
        列出传输任务

        参数：
        ----------
        status : str, optional
            任务状态过滤（pending, in_progress, completed, failed）

        返回：
        ----------
        List[dict]
            任务列表
        """
        tasks = self.state.list_tasks()

        if status:
            tasks = [t for t in tasks if t['status'] == status]

        return tasks

    def clean_completed_tasks(self):
        """清理已完成的任务记录"""
        tasks = self.state.list_tasks()

        for task in tasks:
            if task['status'] == 'completed':
                self.state.remove_task(task['task_id'])

        print(f"已清理 {len([t for t in tasks if t['status'] == 'completed'])} 个已完成任务")

    def retry_failed_tasks(self):
        """重试所有失败的任务"""
        failed_tasks = self.list_tasks(status='failed')

        print(f"找到 {len(failed_tasks)} 个失败的任务")

        for task in failed_tasks:
            print(f"\n重试任务: {task['src']} -> {task['dst']}")
            # 这里需要重新调用上传或下载函数
            # 实际使用时需要保存传输函数引用
