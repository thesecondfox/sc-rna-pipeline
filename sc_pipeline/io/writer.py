"""
数据写入模块

提供数据保存功能
"""

import os
import pandas as pd
import scanpy as sc
from typing import Optional


def save_h5ad(adata, file_path: str, compression: Optional[str] = 'gzip'):
    """
    保存 AnnData 对象为 H5AD 格式

    参数：
    ----------
    adata : AnnData
        要保存的 AnnData 对象
    file_path : str
        输出文件路径
    compression : str, optional
        压缩方法（默认 'gzip'）
    """
    # 确保输出目录存在
    os.makedirs(os.path.dirname(file_path) or '.', exist_ok=True)

    # 保存文件
    adata.write_h5ad(file_path, compression=compression)

    print(f"数据已保存至: {file_path}")


def save_csv(df: pd.DataFrame, file_path: str, **kwargs):
    """
    保存 DataFrame 为 CSV 格式

    参数：
    ----------
    df : pd.DataFrame
        要保存的 DataFrame
    file_path : str
        输出文件路径
    **kwargs
        传递给 pandas.to_csv 的其他参数
    """
    # 确保输出目录存在
    os.makedirs(os.path.dirname(file_path) or '.', exist_ok=True)

    # 保存文件
    df.to_csv(file_path, **kwargs)

    print(f"CSV 已保存至: {file_path}")


def save_figure(fig, file_path: str, dpi: int = 300, **kwargs):
    """
    保存 matplotlib 图形

    参数：
    ----------
    fig : matplotlib.figure.Figure
        要保存的图形对象
    file_path : str
        输出文件路径
    dpi : int
        图像分辨率（默认 300）
    **kwargs
        传递给 savefig 的其他参数
    """
    # 确保输出目录存在
    os.makedirs(os.path.dirname(file_path) or '.', exist_ok=True)

    # 保存图形
    fig.savefig(file_path, dpi=dpi, bbox_inches='tight', **kwargs)

    print(f"图像已保存至: {file_path}")
