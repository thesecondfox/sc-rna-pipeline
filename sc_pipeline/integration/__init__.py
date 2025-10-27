"""
整合模块

提供数据整合、降维和聚类功能
"""

from .integration import data_integration, calculate_heterogeneity, plot_heterogeneity

__all__ = [
    'data_integration',
    'calculate_heterogeneity',
    'plot_heterogeneity'
]
