"""
细胞类型注释模块

提供基于 CellTypist 的细胞类型注释功能
"""

import scanpy as sc
from typing import Optional


def celltypist_annotation(
    adata: sc.AnnData,
    model_path: Optional[str] = None,
    output_path: str = 'annot.h5ad'
):
    """
    使用 CellTypist 进行细胞类型注释

    参数：
    ----------
    adata : sc.AnnData
        输入的 AnnData 对象
    model_path : str
        CellTypist 模型路径（如果为 None，使用默认路径）
    output_path : str
        输出文件路径

    返回：
    ----------
    adata : sc.AnnData
        添加了注释结果的 AnnData 对象
    """
    try:
        import celltypist
        from celltypist import models
    except ImportError:
        raise ImportError(
            "celltypist 未安装。请使用以下命令安装：\n"
            "pip install celltypist"
        )

    if model_path is None:
        model_path = '/share/org/BGI/bgi_zhouy/.celltypist/models/Cells_Intestinal_Tract.pkl'

    print(f"加载 CellTypist 模型: {model_path}")
    model1 = models.Model.load(model=model_path)

    print("开始细胞类型注释...")
    predictions = celltypist.annotate(adata, model=model1, majority_voting=True)

    result_model = predictions.predicted_labels[['predicted_labels', 'majority_voting']].rename(
        columns={'predicted_labels': 'pre1', 'majority_voting': 'pre1_mv'}
    )

    adata.obs = adata.obs.join(result_model, how='left')

    print(f"注释完成，结果保存到 {output_path}")
    adata.write_h5ad(output_path)

    return adata
