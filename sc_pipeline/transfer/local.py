"""
本地文件操作模块

提供本地文件传输和管理功能
"""

import os
import shutil
import hashlib
from pathlib import Path
from typing import Optional, Callable


def calculate_checksum(file_path: str, algorithm: str = 'md5', chunk_size: int = 8192) -> str:
    """
    计算文件校验和

    参数：
    ----------
    file_path : str
        文件路径
    algorithm : str
        哈希算法（md5, sha1, sha256等）
    chunk_size : int
        读取块大小（字节）

    返回：
    ----------
    str
        文件校验和
    """
    hash_obj = hashlib.new(algorithm)

    with open(file_path, 'rb') as f:
        while chunk := f.read(chunk_size):
            hash_obj.update(chunk)

    return hash_obj.hexdigest()


def copy_file(
    src: str,
    dst: str,
    chunk_size: int = 1024 * 1024,
    progress_callback: Optional[Callable] = None,
    verify: bool = True
) -> bool:
    """
    本地文件拷贝（支持进度回调）

    参数：
    ----------
    src : str
        源文件路径
    dst : str
        目标文件路径
    chunk_size : int
        拷贝块大小（字节，默认1MB）
    progress_callback : Callable
        进度回调函数，接收 (已传输字节, 总字节) 参数
    verify : bool
        是否验证校验和（默认True）

    返回：
    ----------
    bool
        是否成功
    """
    # 确保目标目录存在
    os.makedirs(os.path.dirname(dst) or '.', exist_ok=True)

    # 获取文件大小
    file_size = os.path.getsize(src)
    transferred = 0

    # 拷贝文件
    with open(src, 'rb') as src_file:
        with open(dst, 'wb') as dst_file:
            while chunk := src_file.read(chunk_size):
                dst_file.write(chunk)
                transferred += len(chunk)

                # 调用进度回调
                if progress_callback:
                    progress_callback(transferred, file_size)

    # 验证文件完整性
    if verify:
        src_checksum = calculate_checksum(src)
        dst_checksum = calculate_checksum(dst)

        if src_checksum != dst_checksum:
            print(f"警告：校验和不匹配！")
            print(f"  源文件: {src_checksum}")
            print(f"  目标文件: {dst_checksum}")
            return False

    return True


def move_file(src: str, dst: str, verify: bool = True) -> bool:
    """
    移动文件

    参数：
    ----------
    src : str
        源文件路径
    dst : str
        目标文件路径
    verify : bool
        是否验证校验和

    返回：
    ----------
    bool
        是否成功
    """
    # 确保目标目录存在
    os.makedirs(os.path.dirname(dst) or '.', exist_ok=True)

    if verify:
        src_checksum = calculate_checksum(src)

    # 移动文件
    shutil.move(src, dst)

    # 验证
    if verify:
        dst_checksum = calculate_checksum(dst)
        if src_checksum != dst_checksum:
            print(f"警告：移动后校验和不匹配！")
            return False

    return True


def get_file_info(file_path: str) -> dict:
    """
    获取文件信息

    参数：
    ----------
    file_path : str
        文件路径

    返回：
    ----------
    dict
        文件信息字典
    """
    stat = os.stat(file_path)

    return {
        'path': file_path,
        'name': os.path.basename(file_path),
        'size': stat.st_size,
        'size_mb': stat.st_size / 1024 / 1024,
        'modified_time': stat.st_mtime,
        'created_time': stat.st_ctime,
        'is_file': os.path.isfile(file_path),
        'is_dir': os.path.isdir(file_path),
        'exists': os.path.exists(file_path)
    }


def clean_temp_files(directory: str, pattern: str = '*.tmp'):
    """
    清理临时文件

    参数：
    ----------
    directory : str
        目录路径
    pattern : str
        文件匹配模式（默认 *.tmp）
    """
    from glob import glob

    temp_files = glob(os.path.join(directory, pattern))

    for temp_file in temp_files:
        try:
            os.remove(temp_file)
            print(f"已删除临时文件: {temp_file}")
        except Exception as e:
            print(f"删除临时文件失败 {temp_file}: {e}")
