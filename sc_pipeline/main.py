"""
主入口程序

单细胞RNA-seq数据处理流程命令行接口
"""

import os
import argparse
from sc_pipeline.io import read_sample, save_h5ad, save_csv
from sc_pipeline.preprocessing import quality_control
from sc_pipeline.integration import data_integration
from sc_pipeline.annotation import celltypist_annotation
from sc_pipeline.utils import MemoryMonitor, memory_profiler


def get_args():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(description="Single-cell RNA-seq Data Processing Pipeline")

    # 输入输出参数
    parser.add_argument("--sample_info", required=True, help="Path to sample information CSV file")
    parser.add_argument("--output_dir", default="./results", help="Output directory for results")

    # 质控参数
    parser.add_argument("--min_genes", type=int, default=200, help="Minimum genes per cell for QC")
    parser.add_argument("--max_genes", type=int, default=6000, help="Maximum genes per cell for QC")
    parser.add_argument("--max_pct_mito", type=float, default=20.0, help="Maximum percentage of mitochondrial genes for QC")
    parser.add_argument("--doublet_method", default="scrublet", help="Method for doublet detection (scrublet/none)")
    parser.add_argument("--doublet_threshold", type=float, default=0.25, help="Threshold for doublet detection")

    # 整合参数
    parser.add_argument("--integration_method", default="harmony", help="Method for data integration (harmony/combat/none)")
    parser.add_argument("--n_pcs", type=int, default=30, help="Number of principal components")
    parser.add_argument("--n_neighbors", type=int, default=10, help="Number of neighbors for graph construction")
    parser.add_argument("--resolution", type=float, default=1.1, help="Resolution for Leiden clustering")
    parser.add_argument("--gene_num", type=int, default=2000, help="Number of highly variable genes")
    parser.add_argument("--umap_min_dist", type=float, default=0.5, help="Minimum distance for UMAP")

    # 注释参数
    parser.add_argument("--celltypist_model", default=None, help="Path to CellTypist model (optional)")
    parser.add_argument("--run_annotation", action='store_true', help="Run CellTypist annotation")

    # 性能参数
    parser.add_argument("--enable_memory_monitor", action='store_true', help="Enable memory monitoring")
    parser.add_argument("--memory_threshold", type=float, default=80.0, help="Memory usage threshold for warning (%)")

    return parser.parse_args()


@memory_profiler
def main():
    """
    单细胞RNA-seq数据处理主流程
    """
    args = get_args()

    # 启动内存监控（如果启用）
    mem_monitor = None
    if args.enable_memory_monitor:
        print("\n启用内存监控...")
        mem_monitor = MemoryMonitor(
            threshold_percent=args.memory_threshold,
            check_interval=2.0,
            enable_logging=True
        )
        mem_monitor.start_monitoring()

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
    save_h5ad(adata, raw_data_path)

    print(f"\n✓ 读取完成:")
    print(f"  - 细胞数: {adata.n_obs:,}")
    print(f"  - 基因数: {adata.n_vars:,}")
    print(f"  - 数据已保存: {raw_data_path}")

    # 内存检查点
    if mem_monitor:
        mem_monitor.print_summary()

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
    save_h5ad(adata_filtered, filtered_data_path)

    # 保存质控统计表
    qc_stats_path = os.path.join(args.output_dir, "qc_statistics.csv")
    save_csv(qc_stats, qc_stats_path, index=False)

    print(f"\n✓ 质控完成:")
    print(f"  - 过滤后细胞数: {adata_filtered.n_obs:,}")
    print(f"  - .raw 包含基因数: {adata_filtered.raw.n_vars:,} (原始counts)")
    print(f"  - .X 包含基因数: {adata_filtered.n_vars:,} (标准化数据)")
    print(f"  - 质控数据已保存: {filtered_data_path}")
    print(f"  - 质控统计已保存: {qc_stats_path}")

    # 内存检查点
    if mem_monitor:
        mem_monitor.print_summary()

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
    save_h5ad(adata_integrated, integrated_data_path)

    print(f"\n✓ 数据整合完成:")
    print(f"  - 细胞数: {adata_integrated.n_obs:,}")
    print(f"  - .X 中基因数: {adata_integrated.n_vars:,} (高变基因)")
    print(f"  - .raw 中基因数: {adata_integrated.raw.n_vars:,} (完整基因集)")
    print(f"  - 整合数据已保存: {integrated_data_path}")

    # 内存检查点
    if mem_monitor:
        mem_monitor.print_summary()

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

    # 停止内存监控
    if mem_monitor:
        mem_monitor.stop_monitoring()
        print("\n最终内存使用情况：")
        mem_monitor.print_summary()


if __name__ == "__main__":
    main()
