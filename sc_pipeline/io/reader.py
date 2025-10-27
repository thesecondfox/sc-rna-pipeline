"""
数据读取模块

提供单细胞数据读取功能
"""

import os
import scanpy as sc
import pandas as pd
from scipy import io
from scipy.sparse import csr_matrix
import anndata as ad
from tqdm import tqdm


def readbgi(data_dir: str, sample_id: str):
    """
    从 BGI (或 10X) 输出目录读取单细胞矩阵数据，生成 AnnData 对象。

    参数：
    ----------
    data_dir : str
        存放 matrix.mtx(.gz)、barcodes.tsv(.gz)、features.tsv(.gz) 的目录路径
    sample_id : str
        样品 ID，用于为细胞条码添加前缀以避免重复

    返回：
    ----------
    adata : AnnData
        包含表达矩阵与样品标识的 AnnData 对象
    """
    # 自动查找三个文件（支持压缩文件）
    def find_file(name_starts):
        for fn in os.listdir(data_dir):
            if fn.startswith(name_starts) and (fn.endswith(".tsv") or fn.endswith(".tsv.gz") or fn.endswith(".mtx") or fn.endswith(".mtx.gz")):
                return os.path.join(data_dir, fn)
        raise FileNotFoundError(f"未找到 {name_starts} 文件，请检查路径：{data_dir}")

    matrix_path   = find_file("matrix")
    barcodes_path = find_file("barcodes")
    features_path = find_file("features")

    # 读取数据
    matrix = io.mmread(matrix_path)
    barcodes = pd.read_csv(barcodes_path, header=None, sep='\t')
    features = pd.read_csv(features_path, header=None, sep='\t')

    # 给细胞条码加样品前缀
    barcodes[0] = barcodes[0].apply(lambda x: f"{sample_id}_{x}")

    # 转置为细胞×基因矩阵并转换为稀疏格式
    adata = sc.AnnData(X=csr_matrix(matrix.T))
    adata.obs_names = barcodes[0].values
    adata.var_names = features[0].values
    adata.obs["SampleName"] = sample_id  # 只保留 SampleName

    return adata


def read_sample(sample_info_path: str):
    """
    根据样品信息文件批量读取并整合单细胞数据

    参数：
    ----------
    sample_info_path : str
        样品信息 CSV 文件路径
    格式要求：
        必须包含 'Path' 和 'SampleName' 两列，其他列将作为元数据添加到 obs 中。
    返回：
    ----------
    adata : AnnData
        整合后的 AnnData 对象
    """
    # 读取样品信息
    sample_info = pd.read_csv(sample_info_path)

    adata_list = []
    # tqdm 进度条
    for idx, row in tqdm(sample_info.iterrows(), total=len(sample_info), desc="Reading samples"):
        sample_id = row['SampleName']
        data_dir = row['Path']

        # 读取单个样品数据
        adata_sample = readbgi(data_dir, sample_id)

        # 添加元数据
        for col in sample_info.columns:
            if col not in ['Path', 'SampleName']:
                adata_sample.obs[col] = row[col]

        adata_list.append(adata_sample)

    # 合并所有样品
    adata = ad.concat(adata_list, join='outer', merge='same')

    return adata


def read_h5ad(file_path: str):
    """
    读取 H5AD 格式的 AnnData 对象

    参数：
    ----------
    file_path : str
        H5AD 文件路径

    返回：
    ----------
    adata : AnnData
        AnnData 对象
    """
    return sc.read_h5ad(file_path)
