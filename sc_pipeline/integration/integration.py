"""
数据整合模块

提供批次效应校正、降维和聚类功能
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import matplotlib.pyplot as plt
import seaborn as sns
import logging
from scipy.spatial.distance import pdist
from scipy.spatial import ConvexHull
from typing import Optional, List


def data_integration(
    adata: sc.AnnData,
    batch_key: str = "SampleName",
    method: str = "harmony",
    n_pcs: int = 30,
    n_neighbors: int = 10,
    resolution: float = 1.1,
    gene_num: int = 2000,
    umap_min_dist: float = 0.5,
    group_keys: Optional[List[str]] = None,
    output_dir: str = "./integration_results"
) -> sc.AnnData:
    """
    数据整合和可视化

    参数：
    ----------
    adata : sc.AnnData
        质控后的 AnnData 对象
    batch_key : str
        批次效应的列名（默认 "SampleName"）
    method : str
        整合方法，可选 "harmony", "combat", "none"
    n_pcs : int
        主成分数量（默认 30）
    n_neighbors : int
        邻居数量（默认 10）
    resolution : float
        聚类分辨率（默认 1.1）
    gene_num : int
        高变基因数量（默认 2000）
    umap_min_dist : float
        UMAP最小距离参数（默认 0.5）
    group_keys : list
        需要绘制UMAP的分组列名列表（如 ['Study', 'Stage', 'Region']）
    output_dir : str
        输出目录

    返回：
    ----------
    adata : sc.AnnData
        整合后的 AnnData 对象（.raw 包含完整基因集的标准化数据）
    """
    os.makedirs(output_dir, exist_ok=True)

    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    print("=" * 60)
    print("开始数据整合流程...")
    print("=" * 60)

    # 如果未指定分组，使用obs中除了SampleName和技术列外的所有列
    if group_keys is None:
        exclude_cols = [batch_key, 'n_genes_by_counts', 'total_counts', 'pct_counts_mt',
                       'doublet_score', 'is_doublet', 'pass_qc', 'n_genes', 'n_counts']
        group_keys = [col for col in adata.obs.columns if col not in exclude_cols]

    print(f"\n将为以下分组绘制UMAP: {group_keys}")

    # 1. 数据预处理
    print("\n[1/5] 数据标准化和高变基因选择...")

    # 复制adata以避免修改原始数据
    adata = adata.copy()

    # 关键修复：检查数据是否已经标准化
    if adata.X.max() > 100:  # 如果是原始counts
        print("   检测到原始counts，进行标准化...")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    else:
        print("   数据已经标准化")

    # 关键：在选择高变基因之前保存完整的标准化数据到 .raw
    if adata.raw is None:
        adata.raw = adata.copy()
        print(f"   已将完整的标准化数据（{adata.n_vars} 个基因）保存到 .raw")
    else:
        print(f"   .raw 已存在，包含 {adata.raw.n_vars} 个基因")

    # 高变基因选择
    sc.pp.highly_variable_genes(adata, n_top_genes=gene_num, batch_key=batch_key)

    # 重要：使用视图而非 .copy()，这样可以保持与 .raw 的关联
    adata = adata[:, adata.var.highly_variable]

    print(f"   在 .X 中选择了 {adata.n_vars} 个高变基因用于下游分析")
    print(f"   .raw 中保留了完整的 {adata.raw.n_vars} 个基因用于可视化和差异分析")

    # 2. 数据缩放和PCA
    print("\n[2/5] 数据缩放和PCA降维...")

    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=n_pcs)

    print(f"   完成PCA降维 (n_pcs={n_pcs})")

    # 3. 批次效应校正
    if method.lower() == "harmony":
        print(f"\n[3/5] 使用 Harmony 进行批次效应校正...")
        try:
            # 使用 scanpy.external 的 harmony
            sce.pp.harmony_integrate(adata, batch_key, max_iter_harmony=20)
            adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']
            print("   Harmony整合完成")

        except Exception as e:
            print(f"   警告: Harmony 整合失败: {e}")
            print("   将使用未整合的PCA结果")

    elif method.lower() == "combat":
        print(f"\n[3/5] 使用 Combat 进行批次效应校正...")
        sc.pp.combat(adata, key=batch_key)
        print("   Combat整合完成")

    else:
        print(f"\n[3/5] 不进行批次效应校正...")

    # 4. 计算邻居图、UMAP和聚类
    print(f"\n[4/5] 计算邻居图、UMAP和Leiden聚类...")

    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.umap(adata, min_dist=umap_min_dist)
    sc.tl.leiden(adata, resolution=resolution)

    n_clusters = len(adata.obs['leiden'].unique())
    print(f"   UMAP计算完成")
    print(f"   聚类完成，共 {n_clusters} 个cluster")

    # 5. 可视化和异质性分析
    print(f"\n[5/5] 生成可视化图表和异质性分析...")

    # 5.1 绘制不同条件的UMAP
    print("\n   绘制UMAP图...")

    # 设置绘图参数
    sc.set_figure_params(dpi=100, frameon=False, figsize=(8, 6))

    # 绘制leiden聚类
    sc.pl.umap(adata, color='leiden', legend_loc='on data',
               title='Leiden Clustering', show=False, save=False)
    plt.savefig(f"{output_dir}/umap_leiden.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"     ✓ 保存 umap_leiden.png")

    # 绘制批次效应
    sc.pl.umap(adata, color=batch_key, title=f'UMAP by {batch_key}',
               show=False, save=False)
    plt.savefig(f"{output_dir}/umap_{batch_key}.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"     ✓ 保存 umap_{batch_key}.png")

    # 绘制其他分组的UMAP
    for group in group_keys:
        if group in adata.obs.columns:
            try:
                sc.pl.umap(adata, color=group, title=f'UMAP by {group}',
                          show=False, save=False)
                plt.savefig(f"{output_dir}/umap_{group}.png", dpi=300, bbox_inches='tight')
                plt.close()
                print(f"     ✓ 保存 umap_{group}.png")
            except Exception as e:
                print(f"     ✗ 绘制 umap_{group}.png 失败: {e}")

    # 5.2 异质性分析
    heterogeneity_df = calculate_heterogeneity(adata, group_keys)

    # 保存异质性统计表
    if len(heterogeneity_df) > 0:
        heterogeneity_df.to_csv(f"{output_dir}/heterogeneity_stats.csv", index=False)
        print(f"     ✓ 保存 heterogeneity_stats.csv")

    # 5.3 绘制异质性可视化
    if len(heterogeneity_df) > 0:
        print("\n   绘制异质性分析图...")
        plot_heterogeneity(heterogeneity_df, group_keys, output_dir)

    # 打印摘要
    print("\n" + "=" * 60)
    print("数据整合完成！")
    print("=" * 60)
    print(f"整合方法: {method}")
    print(f"细胞总数: {adata.n_obs:,}")
    print(f".X 中基因数（高变基因）: {adata.n_vars:,}")
    print(f".raw 中基因数（完整基因集）: {adata.raw.n_vars:,}")
    print(f"Leiden clusters: {n_clusters}")
    print(f"\n结果保存在: {output_dir}/")
    print("\n生成的文件：")
    print(f"  - umap_leiden.png: Leiden聚类UMAP图")
    print(f"  - umap_{{group}}.png: 各分组条件的UMAP图")
    print(f"  - heterogeneity_stats.csv: 异质性统计表")
    print(f"  - heterogeneity_{{group}}.png: 各分组的异质性分析图")

    return adata


def calculate_heterogeneity(adata: sc.AnnData, group_keys: List[str]) -> pd.DataFrame:
    """
    计算异质性指标

    参数：
    ----------
    adata : sc.AnnData
        包含 UMAP 坐标的 AnnData 对象
    group_keys : list
        需要计算异质性的分组列名列表

    返回：
    ----------
    heterogeneity_df : pd.DataFrame
        异质性统计表
    """
    print("\n   计算异质性指标...")

    heterogeneity_stats = []

    for group in group_keys:
        if group not in adata.obs.columns:
            continue

        categories = adata.obs[group].unique()

        for category in categories:
            mask = adata.obs[group] == category
            if mask.sum() < 10:  # 跳过样本量太小的组
                continue

            # 提取该组的UMAP坐标
            umap_coords = adata[mask].obsm['X_umap']

            # 计算异质性指标
            try:
                # 1. 细胞间平均距离
                distances = pdist(umap_coords, metric='euclidean')
                mean_distance = np.mean(distances)
                std_distance = np.std(distances)

                # 2. 凸包面积（2D UMAP空间的spread）
                hull = ConvexHull(umap_coords)
                hull_area = hull.volume  # 在2D中volume就是面积
            except:
                mean_distance = np.nan
                std_distance = np.nan
                hull_area = np.nan

            # 3. 细胞数量
            n_cells = mask.sum()

            # 4. 占据的cluster数量
            n_clusters_occupied = adata[mask].obs['leiden'].nunique()

            heterogeneity_stats.append({
                'Group': group,
                'Category': category,
                'n_cells': n_cells,
                'mean_distance': mean_distance,
                'std_distance': std_distance,
                'convex_hull_area': hull_area,
                'n_clusters': n_clusters_occupied
            })

    return pd.DataFrame(heterogeneity_stats)


def plot_heterogeneity(heterogeneity_df: pd.DataFrame, group_keys: List[str], output_dir: str):
    """
    绘制异质性分析图

    参数：
    ----------
    heterogeneity_df : pd.DataFrame
        异质性统计表
    group_keys : list
        分组列名列表
    output_dir : str
        输出目录
    """
    for group in group_keys:
        group_data = heterogeneity_df[heterogeneity_df['Group'] == group]

        if len(group_data) == 0:
            continue

        try:
            fig, axes = plt.subplots(2, 2, figsize=(12, 10))
            fig.suptitle(f'Heterogeneity Analysis: {group}', fontsize=16, y=1.02)

            # 细胞数量
            axes[0, 0].barh(group_data['Category'].astype(str), group_data['n_cells'])
            axes[0, 0].set_xlabel('Number of Cells')
            axes[0, 0].set_title('Cell Count')

            # 平均距离
            axes[0, 1].barh(group_data['Category'].astype(str), group_data['mean_distance'])
            axes[0, 1].set_xlabel('Mean Distance')
            axes[0, 1].set_title('Average Cell Distance')

            # 凸包面积
            axes[1, 0].barh(group_data['Category'].astype(str), group_data['convex_hull_area'])
            axes[1, 0].set_xlabel('Convex Hull Area')
            axes[1, 0].set_title('UMAP Space Coverage')

            # cluster数量
            axes[1, 1].barh(group_data['Category'].astype(str), group_data['n_clusters'])
            axes[1, 1].set_xlabel('Number of Clusters')
            axes[1, 1].set_title('Cluster Distribution')

            plt.tight_layout()
            plt.savefig(f"{output_dir}/heterogeneity_{group}.png", dpi=300, bbox_inches='tight')
            plt.close()
            print(f"     ✓ 保存 heterogeneity_{group}.png")
        except Exception as e:
            print(f"     ✗ 绘制 heterogeneity_{group}.png 失败: {e}")
