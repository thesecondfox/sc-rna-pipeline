"""
主入口程序 - 支持断点续传

单细胞RNA-seq数据处理流程命令行接口
"""

import os
import sys
import gc
import argparse
from sc_pipeline.io import read_sample, read_h5ad, save_h5ad, save_csv
from sc_pipeline.preprocessing import quality_control
from sc_pipeline.integration import data_integration
from sc_pipeline.annotation import celltypist_annotation
from sc_pipeline.utils import MemoryMonitor, CheckpointManager


def get_args():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description="单细胞RNA-seq数据处理流程 (支持断点续传)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # 必需参数
    parser.add_argument("--sample_info", required=True,
                       help="样品信息CSV文件路径 (必须包含 Path 和 SampleName 列)")
    parser.add_argument("--output_dir", default="./results",
                       help="输出目录")

    # 断点控制参数
    checkpoint_group = parser.add_argument_group('断点续传控制')
    checkpoint_group.add_argument("--resume", action='store_true',
                                 help="从上次断点继续运行")
    checkpoint_group.add_argument("--reset", action='store_true',
                                 help="重置所有断点，从头开始")
    checkpoint_group.add_argument("--skip-to",
                                 choices=['qc', 'integrate', 'annotate'],
                                 help="跳过到指定步骤（qc/integrate/annotate）")
    checkpoint_group.add_argument("--show-progress", action='store_true',
                                 help="只显示当前进度，不运行")

    # 质控参数
    qc_group = parser.add_argument_group('质控参数')
    qc_group.add_argument("--min_genes", type=int, default=500,
                         help="每个细胞的最低基因数")
    qc_group.add_argument("--max_genes", type=int, default=8000,
                         help="每个细胞的最高基因数")
    qc_group.add_argument("--max_pct_mito", type=float, default=15.0,
                         help="最大线粒体基因百分比")
    qc_group.add_argument("--doublet_method", default="scrublet",
                         choices=["scrublet", "none"],
                         help="双胞检测方法")
    qc_group.add_argument("--doublet_threshold", type=float, default=0.25,
                         help="双胞检测阈值")

    # 整合与聚类参数
    int_group = parser.add_argument_group('整合与聚类参数')
    int_group.add_argument("--integration_method", default="harmony",
                          choices=["harmony", "combat", "none"],
                          help="批次校正方法")
    int_group.add_argument("--n_pcs", type=int, default=50,
                          help="主成分数量")
    int_group.add_argument("--n_neighbors", type=int, default=30,
                          help="邻居数量")
    int_group.add_argument("--resolution", type=float, default=0.8,
                          help="Leiden聚类分辨率")
    int_group.add_argument("--gene_num", type=int, default=3000,
                          help="高变基因数量")
    int_group.add_argument("--umap_min_dist", type=float, default=0.3,
                          help="UMAP最小距离")

    # 注释参数
    ann_group = parser.add_argument_group('细胞类型注释参数')
    ann_group.add_argument("--run_annotation", action='store_true',
                          help="运行CellTypist细胞类型注释")
    ann_group.add_argument("--celltypist_model", default=None,
                          help="CellTypist模型路径")

    return parser.parse_args()


def main():
    """
    单细胞RNA-seq数据处理主流程（支持断点续传）
    """
    args = get_args()

    # 创建输出目录
    os.makedirs(args.output_dir, exist_ok=True)

    # 初始化断点管理器
    ckpt = CheckpointManager(args.output_dir)

    # 处理断点相关命令
    if args.reset:
        ckpt.reset()

    if args.show_progress:
        ckpt.print_status()
        return

    if args.resume or ckpt.checkpoints:
        ckpt.print_status()
        if args.resume:
            print("🔄 从上次断点继续运行...\n")

    # 初始化内存监控
    mem_monitor = MemoryMonitor()
    mem_monitor.checkpoint("流程开始")

    print("\n" + "=" * 100)
    print(" " * 30 + "单细胞RNA-seq数据处理流程 v2.0")
    print("=" * 100)
    print(f"\n📋 配置参数：")
    print(f"  样品信息: {args.sample_info}")
    print(f"  输出目录: {args.output_dir}")
    print(f"  断点文件: {ckpt.checkpoint_file}")
    print(f"\n🔬 质控参数：")
    print(f"  基因数范围: [{args.min_genes}, {args.max_genes}]")
    print(f"  最大线粒体%: {args.max_pct_mito}%")
    print(f"  双胞检测: {args.doublet_method} (阈值={args.doublet_threshold})")
    print(f"\n🧬 分析参数：")
    print(f"  整合方法: {args.integration_method}")
    print(f"  主成分数: {args.n_pcs}")
    print(f"  高变基因数: {args.gene_num}")
    print(f"  聚类分辨率: {args.resolution}")
    print(f"  邻居数: {args.n_neighbors}")
    print("=" * 100 + "\n")

    # ============================================
    # 步骤 1: 读取样品数据
    # ============================================
    raw_data_path = os.path.join(args.output_dir, "01_raw_data.h5ad")

    # 处理跳过命令
    if args.skip_to in ['qc', 'integrate', 'annotate']:
        print("⏭️  跳过步骤1: 数据读取\n")
        ckpt.save_checkpoint('step1_read', file=raw_data_path)

    if not ckpt.is_completed('step1_read'):
        print("\n" + "🔹" * 50)
        print("步骤 1/4: 读取样品数据")
        print("🔹" * 50)

        try:
            adata = read_sample(args.sample_info)

            # 保存原始数据
            save_h5ad(adata, raw_data_path)

            print(f"\n✅ 读取完成:")
            print(f"  细胞数: {adata.n_obs:,}")
            print(f"  基因数: {adata.n_vars:,}")
            print(f"  样品数: {adata.obs['SampleName'].nunique()}")
            print(f"  数据已保存: {raw_data_path}")

            # 保存断点
            ckpt.save_checkpoint('step1_read',
                               file=raw_data_path,
                               n_cells=adata.n_obs,
                               n_genes=adata.n_vars)

            mem_monitor.checkpoint("数据读取完成")

        except Exception as e:
            print(f"\n❌ 读取数据失败: {str(e)}")
            import traceback
            traceback.print_exc()
            sys.exit(1)
    else:
        print("\n✓ 步骤1已完成，加载已有数据...")
        adata = read_h5ad(raw_data_path)
        print(f"  细胞数: {adata.n_obs:,}")
        print(f"  基因数: {adata.n_vars:,}\n")

    # ============================================
    # 步骤 2: 质量控制
    # ============================================
    filtered_data_path = os.path.join(args.output_dir, "02_filtered_data.h5ad")

    # 处理跳过命令
    if args.skip_to in ['integrate', 'annotate']:
        print("⏭️  跳过步骤2: 质量控制\n")
        ckpt.save_checkpoint('step2_qc', file=filtered_data_path)

    if not ckpt.is_completed('step2_qc'):
        print("\n" + "🔹" * 50)
        print("步骤 2/4: 质量控制")
        print("🔹" * 50)

        try:
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
            save_h5ad(adata_filtered, filtered_data_path)

            # 保存质控统计表
            qc_stats_path = os.path.join(args.output_dir, "02_qc_statistics.csv")
            save_csv(qc_stats, qc_stats_path, index=False)

            print(f"\n✅ 质控完成:")
            print(f"  过滤后细胞数: {adata_filtered.n_obs:,}")
            print(f"  数据已保存: {filtered_data_path}")
            print(f"  统计已保存: {qc_stats_path}")

            # 保存断点
            ckpt.save_checkpoint('step2_qc',
                               file=filtered_data_path,
                               n_cells_after=adata_filtered.n_obs)

            mem_monitor.checkpoint("质量控制完成")

            adata = adata_filtered
            del adata_filtered
            gc.collect()

        except Exception as e:
            print(f"\n❌ 质控失败: {str(e)}")
            import traceback
            traceback.print_exc()
            sys.exit(1)
    else:
        print("\n✓ 步骤2已完成，加载已有数据...")
        adata = read_h5ad(filtered_data_path)
        print(f"  过滤后细胞数: {adata.n_obs:,}\n")

    # ============================================
    # 步骤 3: 数据整合
    # ============================================
    integrated_data_path = os.path.join(args.output_dir, "03_integrated_data.h5ad")

    # 处理跳过命令
    if args.skip_to == 'annotate':
        print("⏭️  跳过步骤3: 数据整合\n")
        ckpt.save_checkpoint('step3_integrate', file=integrated_data_path)

    if not ckpt.is_completed('step3_integrate'):
        print("\n" + "🔹" * 50)
        print("步骤 3/4: 数据整合与降维")
        print("🔹" * 50)

        try:
            adata_integrated = data_integration(
                adata,
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
            save_h5ad(adata_integrated, integrated_data_path)

            print(f"\n✅ 数据整合完成:")
            print(f"  细胞数: {adata_integrated.n_obs:,}")
            print(f"  聚类数: {adata_integrated.obs['leiden'].nunique()}")
            print(f"  数据已保存: {integrated_data_path}")

            # 保存断点
            ckpt.save_checkpoint('step3_integrate',
                               file=integrated_data_path,
                               n_clusters=adata_integrated.obs['leiden'].nunique())

            mem_monitor.checkpoint("数据整合完成")

            adata = adata_integrated
            del adata_integrated
            gc.collect()

        except Exception as e:
            print(f"\n❌ 数据整合失败: {str(e)}")
            import traceback
            traceback.print_exc()
            sys.exit(1)
    else:
        print("\n✓ 步骤3已完成，加载已有数据...")
        adata = read_h5ad(integrated_data_path)
        print(f"  整合后细胞数: {adata.n_obs:,}")
        print(f"  聚类数: {adata.obs['leiden'].nunique()}\n")

    # ============================================
    # 步骤 4: 细胞类型注释 (可选)
    # ============================================
    if args.run_annotation:
        annotation_output = os.path.join(args.output_dir, "04_annotated_data.h5ad")

        if not ckpt.is_completed('step4_annotate'):
            print("\n" + "🔹" * 50)
            print("步骤 4/4: 细胞类型注释")
            print("🔹" * 50)

            try:
                adata_annotated = celltypist_annotation(
                    adata,
                    model_path=args.celltypist_model,
                    output_path=annotation_output
                )

                print(f"\n✅ 细胞类型注释完成:")
                print(f"  注释结果已保存: {annotation_output}")

                # 保存断点
                ckpt.save_checkpoint('step4_annotate',
                                   file=annotation_output)

                mem_monitor.checkpoint("细胞类型注释完成")

            except Exception as e:
                print(f"\n⚠️  细胞类型注释失败: {str(e)}")
        else:
            print("\n✓ 步骤4已完成，跳过\n")

    # ============================================
    # 完成总结
    # ============================================
    mem_monitor.checkpoint("流程完成")

    # 保存内存使用报告
    mem_summary = mem_monitor.get_summary()
    mem_report_path = os.path.join(args.output_dir, "memory_usage.csv")
    mem_summary.to_csv(mem_report_path, index=False)

    # 绘制内存使用图
    mem_plot_path = os.path.join(args.output_dir, "memory_usage.png")
    mem_monitor.plot_memory_usage(mem_plot_path)

    print("\n" + "=" * 100)
    print(" " * 40 + "🎉 流程完成！")
    print("=" * 100)
    print(f"\n📁 所有结果已保存至: {args.output_dir}/")
    print("\n📄 生成的文件：")
    print(f"  1. 01_raw_data.h5ad           - 原始数据")
    print(f"  2. 02_filtered_data.h5ad      - 质控后数据")
    print(f"  3. 02_qc_statistics.csv       - 质控统计表")
    print(f"  4. 03_integrated_data.h5ad    - 整合后数据")
    print(f"  5. memory_usage.csv           - 内存使用报告")
    print(f"  6. memory_usage.png           - 内存使用趋势图")

    if args.run_annotation:
        print(f"  7. 04_annotated_data.h5ad     - 细胞类型注释结果")

    print("\n💡 使用提示：")
    print("  • 使用 --resume 参数从断点继续运行")
    print("  • 使用 --reset 参数重置所有断点")
    print("  • 使用 --show-progress 查看当前进度")
    print("  • 使用 --skip-to [qc|integrate|annotate] 跳过到指定步骤")

    print("\n📊 流程统计：")
    print(f"  峰值内存: {mem_monitor.get_memory_usage():.2f} GB")
    print(f"  完成进度: {ckpt.get_progress_percentage():.0f}%")

    print("=" * 100 + "\n")


if __name__ == "__main__":
    main()
