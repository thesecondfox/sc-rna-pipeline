"""
质量控制模块

提供单细胞数据质控功能，包括：
- QC 指标计算
- 双胞检测
- 细胞过滤
"""

import scanpy as sc
import pandas as pd
from typing import Union, Tuple
from tqdm import tqdm
import scrublet as scr


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

    # 过滤数据并正确保留 .raw
    adata_filtered = adata[adata.obs['pass_qc']].copy()

    # 关键修复：在标准化前先保存原始counts到 .raw
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
