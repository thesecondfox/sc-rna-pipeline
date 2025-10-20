import os
import scanpy as sc
import pandas as pd
from scipy import io
from scipy.sparse import csr_matrix
import anndata as ad
import numpy as np
from typing import Union, Tuple
from tqdm import tqdm
import scrublet as scr


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
    from tqdm import tqdm
    
    # 读取样品信息
    sample_info = pd.read_csv(sample_info_path)
    
    adata_list = []
    # ✅ tqdm 进度条
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


def quality_control(
    adata: sc.AnnData,
    min_genes: int = 200,
    max_genes: int = 6000,
    max_pct_mito: float = 20,
    doublet_method: str = "scrublet",
    doublet_threshold: float = 0.25,
    return_stats: bool = True
) -> Union[sc.AnnData, Tuple[sc.AnnData, pd.DataFrame]]:
    """
    单细胞数据质控流程
    
    参数：
    ----------
    adata : sc.AnnData
        输入的 AnnData 对象
    min_genes : int
        最低基因数阈值（默认 200）
    max_genes : int
        最高基因数阈值（默认 6000）
    max_pct_mito : float
        最大线粒体比例阈值（默认 20%）
    doublet_method : str
        双胞去除方法，可选 "scrublet" 或 "none"
    doublet_threshold : float
        双胞预测阈值（默认 0.25）
    return_stats : bool
        是否返回质控统计表（默认 True）
    
    返回：
    ----------
    adata_filtered : sc.AnnData
        质控后的 AnnData 对象（.raw 包含原始counts）
    qc_stats : pd.DataFrame (可选)
        每个样品的质控统计表
    """
    from tqdm import tqdm
    import scrublet as scr
    
    print("=" * 60)
    print("开始质控流程...")
    print("=" * 60)
    
    # 1. 计算QC指标
    print("\n[1/4] 计算质控指标...")
    
    # 识别线粒体基因（支持人和鼠）
    adata.var['mt'] = adata.var_names.str.startswith('MT-') | adata.var_names.str.startswith('mt-')
    
    # 计算QC metrics
    sc.pp.calculate_qc_metrics(
        adata, 
        qc_vars=['mt'], 
        percent_top=None, 
        log1p=False, 
        inplace=True
    )
    
    # 记录过滤前的细胞数
    samples = adata.obs['SampleName'].unique()
    qc_stats_list = []
    
    for sample in samples:
        sample_mask = adata.obs['SampleName'] == sample
        n_cells_before = sample_mask.sum()
        qc_stats_list.append({
            'SampleName': sample,
            'cells_before_qc': n_cells_before
        })
    
    # 2. 去除双胞
    if doublet_method.lower() == "scrublet":
        print(f"\n[2/4] 使用 Scrublet 检测双胞 (阈值={doublet_threshold})...")
        
        doublet_scores = []
        doublet_labels = []
        
        for sample in tqdm(samples, desc="Detecting doublets"):
            sample_mask = adata.obs['SampleName'] == sample
            adata_sample = adata[sample_mask].copy()
            
            # 运行 Scrublet
            scrub = scr.Scrublet(adata_sample.X)
            doublet_score, predicted_doublet = scrub.scrub_doublets(
                min_counts=2,
                min_cells=3,
                min_gene_variability_pctl=85,
                n_prin_comps=30
            )
            
            # 使用自定义阈值
            predicted_doublet = doublet_score > doublet_threshold
            
            doublet_scores.extend(doublet_score)
            doublet_labels.extend(predicted_doublet)
        
        adata.obs['doublet_score'] = doublet_scores
        adata.obs['is_doublet'] = doublet_labels
        
        n_doublets = sum(doublet_labels)
        print(f"   检测到 {n_doublets} 个双胞 ({n_doublets/len(adata)*100:.2f}%)")
    else:
        print("\n[2/4] 跳过双胞检测...")
        adata.obs['is_doublet'] = False
    
    # 3. 应用过滤条件
    print(f"\n[3/4] 应用过滤条件...")
    print(f"   - 最低基因数: {min_genes}")
    print(f"   - 最高基因数: {max_genes}")
    print(f"   - 最大线粒体比例: {max_pct_mito}%")
    
    # 创建过滤mask
    adata.obs['pass_qc'] = (
        (adata.obs['n_genes_by_counts'] >= min_genes) &
        (adata.obs['n_genes_by_counts'] <= max_genes) &
        (adata.obs['pct_counts_mt'] <= max_pct_mito) &
        (~adata.obs['is_doublet'])
    )
    
    # 4. 统计每个样品的过滤情况
    print("\n[4/4] 统计过滤结果...")
    
    for i, sample in enumerate(samples):
        sample_mask = adata.obs['SampleName'] == sample
        sample_data = adata.obs[sample_mask]
        
        n_cells_before = qc_stats_list[i]['cells_before_qc']
        n_cells_after = sample_data['pass_qc'].sum()
        n_filtered = n_cells_before - n_cells_after
        pct_filtered = n_filtered / n_cells_before * 100
        
        # 统计各种过滤原因
        n_low_genes = (sample_data['n_genes_by_counts'] < min_genes).sum()
        n_high_genes = (sample_data['n_genes_by_counts'] > max_genes).sum()
        n_high_mito = (sample_data['pct_counts_mt'] > max_pct_mito).sum()
        n_doublet = sample_data['is_doublet'].sum()
        
        qc_stats_list[i].update({
            'cells_after_qc': n_cells_after,
            'cells_filtered': n_filtered,
            'pct_filtered': pct_filtered,
            'n_low_genes': n_low_genes,
            'n_high_genes': n_high_genes,
            'n_high_mito': n_high_mito,
            'n_doublet': n_doublet,
            'mean_genes': sample_data[sample_data['pass_qc']]['n_genes_by_counts'].mean(),
            'median_genes': sample_data[sample_data['pass_qc']]['n_genes_by_counts'].median(),
            'mean_counts': sample_data[sample_data['pass_qc']]['total_counts'].mean(),
            'mean_pct_mito': sample_data[sample_data['pass_qc']]['pct_counts_mt'].mean()
        })
    
    qc_stats_df = pd.DataFrame(qc_stats_list)
    
    # ✅ 过滤数据并正确保留 .raw
    adata_filtered = adata[adata.obs['pass_qc']].copy()
    
    # ✅ 关键修复：在标准化前先保存原始counts到 .raw
    if adata_filtered.X.max() > 100:  # 如果是原始counts
        # 先保存原始counts到 .raw
        adata_filtered.raw = adata_filtered.copy()
        print(f"   已将原始counts（{adata_filtered.n_vars} 个基因）保存到 .raw")
        
        # 再对 .X 进行标准化
        sc.pp.normalize_total(adata_filtered, target_sum=1e4)
        sc.pp.log1p(adata_filtered)
        print("   .X 已标准化（log-normalized）")
    else:
        # 如果数据已经标准化
        if adata_filtered.raw is None:
            adata_filtered.raw = adata_filtered.copy()
            print("   数据已标准化，创建 .raw 副本")
        else:
            # 如果原始数据有 .raw，需要同步过滤
            adata_filtered.raw = adata.raw[adata.obs['pass_qc']].copy()
            print("   已同步过滤 .raw 数据")
    
    # 打印总体统计
    print("\n" + "=" * 60)
    print("质控完成！")
    print("=" * 60)
    print(f"过滤前总细胞数: {len(adata):,}")
    print(f"过滤后总细胞数: {len(adata_filtered):,}")
    print(f"总过滤比例: {(len(adata) - len(adata_filtered))/len(adata)*100:.2f}%")
    print(f".raw 包含: {adata_filtered.raw.n_vars} 个基因（原始counts）")
    print(f".X 包含: {adata_filtered.n_vars} 个基因（标准化数据）")
    print("\n各样品质控统计：")
    print(qc_stats_df[['SampleName', 'cells_before_qc', 'cells_after_qc', 'pct_filtered']].to_string(index=False))
    
    # 标记异常样品（过滤比例>50%）
    abnormal_samples = qc_stats_df[qc_stats_df['pct_filtered'] > 50]
    if len(abnormal_samples) > 0:
        print("\n⚠️  警告：以下样品过滤比例超过50%，建议检查：")
        print(abnormal_samples[['SampleName', 'cells_before_qc', 'cells_after_qc', 'pct_filtered']].to_string(index=False))
    
    if return_stats:
        return adata_filtered, qc_stats_df
    else:
        return adata_filtered


def data_integration(
    adata: sc.AnnData,
    batch_key: str = "SampleName",
    method: str = "harmony",
    n_pcs: int = 30,
    n_neighbors: int = 10,
    resolution: float = 1.1,
    gene_num: int = 2000,
    umap_min_dist: float = 0.5,
    group_keys: list = None,
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
    import os
    import matplotlib.pyplot as plt
    import seaborn as sns
    import scanpy.external as sce
    import logging
    from scipy.spatial.distance import pdist
    from scipy.spatial import ConvexHull
    
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
    
    # ✅ 关键修复：检查数据是否已经标准化
    if adata.X.max() > 100:  # 如果是原始counts
        print("   检测到原始counts，进行标准化...")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    else:
        print("   数据已经标准化")
    
    # ✅ 关键：在选择高变基因之前保存完整的标准化数据到 .raw
    if adata.raw is None:
        adata.raw = adata.copy()
        print(f"   已将完整的标准化数据（{adata.n_vars} 个基因）保存到 .raw")
    else:
        print(f"   .raw 已存在，包含 {adata.raw.n_vars} 个基因")
    
    # 高变基因选择
    sc.pp.highly_variable_genes(adata, n_top_genes=gene_num, batch_key=batch_key)
    
    # ✅ 重要：使用视图而非 .copy()，这样可以保持与 .raw 的关联
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
    
    heterogeneity_df = pd.DataFrame(heterogeneity_stats)
    
    # 保存异质性统计表
    if len(heterogeneity_df) > 0:
        heterogeneity_df.to_csv(f"{output_dir}/heterogeneity_stats.csv", index=False)
        print(f"     ✓ 保存 heterogeneity_stats.csv")
    
    # 5.3 绘制异质性可视化
    if len(heterogeneity_df) > 0:
        print("\n   绘制异质性分析图...")
        
        for group in group_keys:
            if group not in adata.obs.columns:
                continue
            
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


def celltypist_annotation(adata: sc.AnnData, model_path: str = None, output_path: str = 'annot.h5ad'):
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
    import celltypist
    from celltypist import models
    
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


def get_args():
    import argparse

    parser = argparse.ArgumentParser(description="Single-cell RNA-seq Data Processing Pipeline")
    parser.add_argument("--sample_info", required=True, help="Path to sample information CSV file")
    parser.add_argument("--output_dir", default="./results", help="Output directory for results")
    parser.add_argument("--min_genes", type=int, default=200, help="Minimum genes per cell for QC")
    parser.add_argument("--max_genes", type=int, default=6000, help="Maximum genes per cell for QC")
    parser.add_argument("--max_pct_mito", type=float, default=20.0, help="Maximum percentage of mitochondrial genes for QC")
    parser.add_argument("--doublet_method", default="scrublet", help="Method for doublet detection (scrublet/none)")
    parser.add_argument("--doublet_threshold", type=float, default=0.25, help="Threshold for doublet detection")
    parser.add_argument("--integration_method", default="harmony", help="Method for data integration (harmony/combat/none)")
    parser.add_argument("--n_pcs", type=int, default=30, help="Number of principal components")
    parser.add_argument("--n_neighbors", type=int, default=10, help="Number of neighbors for graph construction")
    parser.add_argument("--resolution", type=float, default=1.1, help="Resolution for Leiden clustering")
    parser.add_argument("--gene_num", type=int, default=2000, help="Number of highly variable genes")
    parser.add_argument("--umap_min_dist", type=float, default=0.5, help="Minimum distance for UMAP")
    parser.add_argument("--celltypist_model", default=None, help="Path to CellTypist model (optional)")
    parser.add_argument("--run_annotation", action='store_true', help="Run CellTypist annotation")
    
    return parser.parse_args()


def main():
    """
    单细胞RNA-seq数据处理主流程
    """
    args = get_args()
    
    # 创建输出目录
    os.makedirs(args.output_dir, exist_ok=True)
    
    print("\n" + "=" * 80)
    print("单细胞RNA-seq数据处理流程")
    print("=" * 80)
    print(f"\n配置参数：")
    print(f"  - 样品信息文件: {args.sample_info}")
    print(f"  - 输出目录: {args.output_dir}")
    print(f"  - 质控参数: min_genes={args.min_genes}, max_genes={args.max_genes}, max_pct_mito={args.max_pct_mito}")
    print(f"  - 双胞检测: method={args.doublet_method}, threshold={args.doublet_threshold}")
    print(f"  - 整合方法: {args.integration_method}")
    print(f"  - 聚类参数: n_pcs={args.n_pcs}, resolution={args.resolution}")
    print("=" * 80 + "\n")
    
    # ============================================
    # 步骤 1: 读取样品数据
    # ============================================
    print("\n" + ">" * 80)
    print("步骤 1: 读取样品数据")
    print(">" * 80)
    
    adata = read_sample(args.sample_info)
    
    # 保存原始数据
    raw_data_path = os.path.join(args.output_dir, "raw_data.h5ad")
    adata.write(raw_data_path)
    
    print(f"\n✓ 读取完成:")
    print(f"  - 细胞数: {adata.n_obs:,}")
    print(f"  - 基因数: {adata.n_vars:,}")
    print(f"  - 数据已保存: {raw_data_path}")
    
    # ============================================
    # 步骤 2: 质量控制
    # ============================================
    print("\n" + ">" * 80)
    print("步骤 2: 质量控制")
    print(">" * 80)
    
    adata_filtered, qc_stats = quality_control(
        adata,
        min_genes=args.min_genes,
        max_genes=args.max_genes,
        max_pct_mito=args.max_pct_mito,
        doublet_method=args.doublet_method,
        doublet_threshold=args.doublet_threshold,
        return_stats=True
    )
    
    # 保存质控后的数据
    filtered_data_path = os.path.join(args.output_dir, "filtered_data.h5ad")
    adata_filtered.write(filtered_data_path)
    
    # 保存质控统计表
    qc_stats_path = os.path.join(args.output_dir, "qc_statistics.csv")
    qc_stats.to_csv(qc_stats_path, index=False)
    
    print(f"\n✓ 质控完成:")
    print(f"  - 过滤后细胞数: {adata_filtered.n_obs:,}")
    print(f"  - .raw 包含基因数: {adata_filtered.raw.n_vars:,} (原始counts)")
    print(f"  - .X 包含基因数: {adata_filtered.n_vars:,} (标准化数据)")
    print(f"  - 质控数据已保存: {filtered_data_path}")
    print(f"  - 质控统计已保存: {qc_stats_path}")
    
    # ============================================
    # 步骤 3: 数据整合
    # ============================================
    print("\n" + ">" * 80)
    print("步骤 3: 数据整合与降维")
    print(">" * 80)
    
    adata_integrated = data_integration(
        adata_filtered,
        batch_key="SampleName",
        method=args.integration_method,
        n_pcs=args.n_pcs,
        n_neighbors=args.n_neighbors,
        resolution=args.resolution,
        gene_num=args.gene_num,
        umap_min_dist=args.umap_min_dist,
        output_dir=args.output_dir
    )
    
    # 保存整合后的数据
    integrated_data_path = os.path.join(args.output_dir, "integrated_data.h5ad")
    adata_integrated.write(integrated_data_path)
    
    print(f"\n✓ 数据整合完成:")
    print(f"  - 细胞数: {adata_integrated.n_obs:,}")
    print(f"  - .X 中基因数: {adata_integrated.n_vars:,} (高变基因)")
    print(f"  - .raw 中基因数: {adata_integrated.raw.n_vars:,} (完整基因集)")
    print(f"  - 整合数据已保存: {integrated_data_path}")
    
    # ============================================
    # 步骤 4: 细胞类型注释 (可选)
    # ============================================
    if args.run_annotation:
        print("\n" + ">" * 80)
        print("步骤 4: 细胞类型注释")
        print(">" * 80)
        
        annotation_output = os.path.join(args.output_dir, "annotated_data.h5ad")
        
        try:
            adata_annotated = celltypist_annotation(
                adata_integrated,
                model_path=args.celltypist_model,
                output_path=annotation_output
            )
            
            print(f"\n✓ 细胞类型注释完成:")
            print(f"  - 注释结果已保存: {annotation_output}")
            
        except Exception as e:
            print(f"\n✗ 细胞类型注释失败: {e}")
            print("  请检查 CellTypist 模型路径是否正确")
    
    # ============================================
    # 完成总结
    # ============================================
    print("\n" + "=" * 80)
    print("流程完成！")
    print("=" * 80)
    print(f"\n所有结果已保存至: {args.output_dir}/")
    print("\n生成的文件：")
    print(f"  1. raw_data.h5ad - 原始数据")
    print(f"  2. filtered_data.h5ad - 质控后数据 (.raw 包含原始counts)")
    print(f"  3. qc_statistics.csv - 质控统计表")
    print(f"  4. integrated_data.h5ad - 整合后数据 (.raw 包含完整基因集)")
    print(f"  5. umap_*.png - 各种UMAP可视化图")
    print(f"  6. heterogeneity_*.png - 异质性分析图")
    print(f"  7. heterogeneity_stats.csv - 异质性统计表")
    
    if args.run_annotation:
        print(f"  8. annotated_data.h5ad - 细胞类型注释结果")
    
    print("\n提示：")
    print("  - 使用 scanpy.read_h5ad() 可以读取任何 .h5ad 文件")
    print("  - adata.raw 包含完整基因集，可用于差异表达分析")
    print("  - adata.X 包含高变基因，已用于降维和聚类")
    print("=" * 80 + "\n")


if __name__ == "__main__":
    main()


"""
使用示例：

基础运行（不做细胞类型注释）：
python sc_pipeline.py \
    --sample_info samples.csv \
    --output_dir ./results \
    --min_genes 200 \
    --max_genes 6000 \
    --max_pct_mito 20 \
    --doublet_method scrublet \
    --doublet_threshold 0.25 \
    --integration_method harmony \
    --n_pcs 30 \
    --resolution 1.1 \
    --gene_num 2000

包含细胞类型注释：
python sc_pipeline.py \
    --sample_info samples.csv \
    --output_dir ./results \
    --min_genes 200 \
    --max_genes 6000 \
    --max_pct_mito 20 \
    --doublet_method scrublet \
    --doublet_threshold 0.25 \
    --integration_method harmony \
    --run_annotation \
    --celltypist_model /path/to/your/model.pkl

参数说明：
- sample_info: 包含 'Path' 和 'SampleName' 列的 CSV 文件
- output_dir: 输出目录
- min_genes/max_genes: 质控的基因数阈值
- max_pct_mito: 线粒体基因比例阈值
- doublet_method: 双胞检测方法 (scrublet/none)
- doublet_threshold: 双胞检测阈值
- integration_method: 批次校正方法 (harmony/combat/none)
- n_pcs: 主成分数量
- resolution: Leiden 聚类分辨率
- gene_num: 高变基因数量
- run_annotation: 是否运行细胞类型注释
- celltypist_model: CellTypist 模型路径
"""
