"""
单细胞RNA-seq数据分析流程

模块化的单细胞RNA测序数据处理工具包，支持断点续传和内存监控

主要功能：
- 数据读取和批量处理
- 质量控制和过滤
- 批次效应校正
- 降维和聚类
- 细胞类型注释
- 断点续传（流程中断后可继续）
- 内存监控
"""

__version__ = '2.0.0'
__author__ = 'BGI'

# IO模块
from .io import (
    readbgi,
    read_sample,
    read_h5ad,
    save_h5ad,
    save_csv,
    save_figure
)

# 预处理模块
from .preprocessing import quality_control

# 整合模块
from .integration import (
    data_integration,
    calculate_heterogeneity,
    plot_heterogeneity
)

# 注释模块
from .annotation import celltypist_annotation

# 工具模块
from .utils import (
    MemoryMonitor,
    log_memory,
    memory_profiler,
    CheckpointManager
)

__all__ = [
    # Version
    '__version__',
    '__author__',

    # IO
    'readbgi',
    'read_sample',
    'read_h5ad',
    'save_h5ad',
    'save_csv',
    'save_figure',

    # Preprocessing
    'quality_control',

    # Integration
    'data_integration',
    'calculate_heterogeneity',
    'plot_heterogeneity',

    # Annotation
    'celltypist_annotation',

    # Utils
    'MemoryMonitor',
    'log_memory',
    'memory_profiler',
    'CheckpointManager'
]
