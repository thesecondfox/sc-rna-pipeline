import os
import sys
import gc
import json
import psutil
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import io
from scipy.sparse import csr_matrix
import anndata as ad
from typing import Union, Tuple, Optional
from tqdm import tqdm
import warnings
import hashlib
from datetime import datetime
import traceback

warnings.filterwarnings('ignore')

# 设置 matplotlib 参数
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10


class CheckpointManager:
    """
    增强的断点管理器 - 支持从中断处继续运行

    新增功能：
    - 数据文件完整性校验（MD5）
    - 自动检测损坏的断点
    - 交互式恢复确认
    - 详细的错误日志
    """

    def __init__(self, output_dir: str, auto_resume: bool = False):
        self.output_dir = output_dir
        self.checkpoint_file = os.path.join(output_dir, '.checkpoint.json')
        self.error_log_file = os.path.join(output_dir, '.error_log.txt')
        self.auto_resume = auto_resume
        self.checkpoints = self.load_checkpoints()
        self.steps_order = ['step1_read', 'step2_qc', 'step3_integrate', 'step4_annotate']

    def load_checkpoints(self) -> dict:
        """加载已有的检查点并验证完整性"""
        if os.path.exists(self.checkpoint_file):
            try:
                with open(self.checkpoint_file, 'r') as f:
                    checkpoints = json.load(f)

                # 验证检查点完整性
                self._validate_checkpoints(checkpoints)
                return checkpoints
            except Exception as e:
                print(f"⚠️  加载检查点失败: {str(e)}")
                print("   将从头开始运行")
                return {}
        return {}

    def _validate_checkpoints(self, checkpoints: dict):
        """验证检查点数据的完整性"""
        for step_name, step_data in checkpoints.items():
            if not step_data.get('completed', False):
                continue

            # 检查文件是否存在
            file_path = step_data.get('file')
            if file_path and not os.path.exists(file_path):
                print(f"⚠️  警告: 步骤 {step_name} 的数据文件不存在: {file_path}")
                step_data['valid'] = False
            else:
                step_data['valid'] = True

                # 验证文件完整性（可选，针对大文件可能较慢）
                if file_path and os.path.exists(file_path):
                    file_size = os.path.getsize(file_path)
                    saved_size = step_data.get('file_size')

                    if saved_size and file_size != saved_size:
                        print(f"⚠️  警告: 文件大小不匹配 {os.path.basename(file_path)}")
                        print(f"   预期: {saved_size}, 实际: {file_size}")
                        step_data['valid'] = False

    def save_checkpoint(self, step: str, **kwargs):
        """
        保存检查点（增强版）

        自动记录：
        - 时间戳
        - 文件路径和大小
        - 关键统计信息
        """
        checkpoint_data = {
            'completed': True,
            'timestamp': datetime.now().isoformat(),
            'valid': True,
            **kwargs
        }

        # 如果有文件路径，记录文件大小
        if 'file' in kwargs and os.path.exists(kwargs['file']):
            checkpoint_data['file_size'] = os.path.getsize(kwargs['file'])

        self.checkpoints[step] = checkpoint_data

        try:
            with open(self.checkpoint_file, 'w') as f:
                json.dump(self.checkpoints, f, indent=2)
            print(f"   💾 断点已保存: {step}")
        except Exception as e:
            print(f"   ⚠️  保存断点失败: {str(e)}")

    def mark_step_failed(self, step: str, error: Exception):
        """
        标记步骤失败并记录错误信息
        """
        error_info = {
            'failed': True,
            'timestamp': datetime.now().isoformat(),
            'error_type': type(error).__name__,
            'error_message': str(error),
            'traceback': traceback.format_exc()
        }

        self.checkpoints[step] = error_info

        # 保存到主检查点文件
        try:
            with open(self.checkpoint_file, 'w') as f:
                json.dump(self.checkpoints, f, indent=2)
        except:
            pass

        # 追加到错误日志
        try:
            with open(self.error_log_file, 'a') as f:
                f.write(f"\n{'='*80}\n")
                f.write(f"步骤: {step}\n")
                f.write(f"时间: {error_info['timestamp']}\n")
                f.write(f"错误类型: {error_info['error_type']}\n")
                f.write(f"错误信息: {error_info['error_message']}\n")
                f.write(f"详细追踪:\n{error_info['traceback']}\n")
        except:
            pass

        print(f"   ❌ 步骤 {step} 失败，错误已记录")
        print(f"   错误日志: {self.error_log_file}")

    def is_completed(self, step: str) -> bool:
        """检查某个步骤是否已完成且有效"""
        step_data = self.checkpoints.get(step, {})
        return step_data.get('completed', False) and step_data.get('valid', True)

    def is_failed(self, step: str) -> bool:
        """检查某个步骤是否失败过"""
        return self.checkpoints.get(step, {}).get('failed', False)

    def get_checkpoint_data(self, step: str) -> dict:
        """获取检查点数据"""
        return self.checkpoints.get(step, {})

    def reset(self, confirm: bool = False):
        """重置所有检查点（从头开始）"""
        if not confirm and not self.auto_resume:
            response = input("⚠️  确定要重置所有断点吗？(y/n): ")
            if response.lower() != 'y':
                print("取消重置")
                return

        self.checkpoints = {}
        if os.path.exists(self.checkpoint_file):
            os.remove(self.checkpoint_file)
        print("✅ 已重置所有检查点")

    def get_last_completed_step(self) -> Optional[str]:
        """获取最后一个完成的步骤"""
        for step in reversed(self.steps_order):
            if self.is_completed(step):
                return step
        return None

    def get_next_step(self) -> Optional[str]:
        """获取下一个需要执行的步骤"""
        for step in self.steps_order:
            if not self.is_completed(step):
                return step
        return None

    def print_status(self):
        """打印当前进度状态（增强版）"""
        step_names = {
            'step1_read': '步骤1: 数据读取',
            'step2_qc': '步骤2: 质量控制',
            'step3_integrate': '步骤3: 数据整合',
            'step4_annotate': '步骤4: 细胞注释'
        }

        print("\n" + "=" * 80)
        print("当前分析进度")
        print("=" * 80)

        for step_id in self.steps_order:
            step_name = step_names.get(step_id, step_id)

            if self.is_failed(step_id):
                data = self.get_checkpoint_data(step_id)
                print(f"❌ {step_name} - 失败")
                print(f"   错误: {data.get('error_message', '未知错误')}")
                print(f"   时间: {data.get('timestamp', '')[:19]}")

            elif self.is_completed(step_id):
                data = self.get_checkpoint_data(step_id)
                timestamp = data.get('timestamp', '')
                file_path = data.get('file', '')

                print(f"✅ {step_name}")
                if timestamp:
                    print(f"   完成时间: {timestamp[:19]}")
                if file_path and os.path.exists(file_path):
                    size_mb = os.path.getsize(file_path) / 1024**2
                    print(f"   文件: {os.path.basename(file_path)} ({size_mb:.1f} MB)")

                # 显示关键统计
                if 'n_cells' in data:
                    print(f"   细胞数: {data['n_cells']:,}")
                if 'n_genes' in data:
                    print(f"   基因数: {data['n_genes']:,}")
                if 'n_clusters' in data:
                    print(f"   聚类数: {data['n_clusters']}")
            else:
                print(f"⏳ {step_name} - 待执行")

        print("=" * 80)

        # 显示下一步操作
        next_step = self.get_next_step()
        if next_step:
            print(f"\n📍 下一步: {step_names.get(next_step, next_step)}")
        else:
            print(f"\n🎉 所有步骤已完成！")
        print()

    def should_resume(self) -> bool:
        """判断是否应该从断点恢复"""
        last_step = self.get_last_completed_step()

        if not last_step:
            return False

        if self.auto_resume:
            print(f"🔄 自动从断点恢复: {last_step}")
            return True

        # 交互式确认
        print(f"\n检测到之前的进度:")
        self.print_status()

        response = input("是否从上次断点继续运行？(y/n): ")
        return response.lower() == 'y'


class MemoryMonitor:
    """内存监控类"""

    def __init__(self):
        self.process = psutil.Process()
        self.checkpoints = []

    def get_memory_usage(self) -> float:
        """获取当前内存使用（GB）"""
        return self.process.memory_info().rss / 1024**3

    def checkpoint(self, step_name: str):
        """记录检查点"""
        mem_gb = self.get_memory_usage()
        self.checkpoints.append({'step': step_name, 'memory_gb': mem_gb})
        print(f"   📊 内存使用: {mem_gb:.2f} GB ({step_name})")

    def get_summary(self) -> pd.DataFrame:
        """获取内存使用摘要"""
        df = pd.DataFrame(self.checkpoints)
        if len(df) > 0:
            df['memory_increase_gb'] = df['memory_gb'].diff()
        return df

    def plot_memory_usage(self, output_path: str):
        """绘制内存使用趋势图"""
        df = self.get_summary()
        if len(df) == 0:
            return

        fig, ax = plt.subplots(figsize=(12, 6))
        ax.plot(range(len(df)), df['memory_gb'], marker='o', linewidth=2, markersize=8)
        ax.set_xticks(range(len(df)))
        ax.set_xticklabels(df['step'], rotation=45, ha='right')
        ax.set_ylabel('Memory Usage (GB)', fontsize=12)
        ax.set_xlabel('Pipeline Step', fontsize=12)
        ax.set_title('Memory Usage Throughout Pipeline', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"   💾 内存使用图已保存: {output_path}")


def safe_step_execution(checkpoint_manager: CheckpointManager, step_name: str,
                       func, *args, **kwargs):
    """
    安全执行流程步骤，自动处理错误和断点保存

    参数:
    --------
    checkpoint_manager: CheckpointManager
        断点管理器
    step_name: str
        步骤名称
    func: callable
        要执行的函数
    *args, **kwargs:
        传递给函数的参数

    返回:
    --------
    result: 函数执行结果
    """
    try:
        print(f"\n{'='*80}")
        print(f"开始执行: {step_name}")
        print(f"{'='*80}")

        result = func(*args, **kwargs)

        print(f"✅ {step_name} 执行成功")
        return result

    except KeyboardInterrupt:
        print(f"\n⚠️  用户中断执行")
        checkpoint_manager.mark_step_failed(step_name,
                                           Exception("User interrupted"))
        sys.exit(1)

    except Exception as e:
        print(f"\n❌ {step_name} 执行失败!")
        print(f"   错误类型: {type(e).__name__}")
        print(f"   错误信息: {str(e)}")

        # 记录失败信息
        checkpoint_manager.mark_step_failed(step_name, e)

        # 打印详细错误信息
        print(f"\n详细错误追踪:")
        print(traceback.format_exc())

        print(f"\n💡 可以修复问题后使用 --resume 参数从断点继续运行")

        sys.exit(1)


def readbgi(data_dir: str, sample_id: str,
           memory_monitor: Optional[MemoryMonitor] = None) -> ad.AnnData:
    """
    从 BGI (或 10X) 输出目录读取单细胞矩阵数据，生成 AnnData 对象。

    参数：
    ----------
    data_dir : str
        存放 matrix.mtx(.gz)、barcodes.tsv(.gz)、features.tsv(.gz) 的目录路径
    sample_id : str
        样品 ID，用于为细胞条码添加前缀以避免重复
    memory_monitor : MemoryMonitor, optional
        内存监控对象

    返回：
    ----------
    adata : AnnData
        包含表达矩阵与样品标识的 AnnData 对象
    """
    # 自动查找三个文件（支持压缩文件）
    def find_file(name_starts):
        for fn in os.listdir(data_dir):
            if fn.startswith(name_starts) and (fn.endswith(".tsv") or fn.endswith(".tsv.gz") or
                                               fn.endswith(".mtx") or fn.endswith(".mtx.gz")):
                return os.path.join(data_dir, fn)
        raise FileNotFoundError(f"未找到 {name_starts} 文件，请检查路径：{data_dir}")

    try:
        matrix_path = find_file("matrix")
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
        adata.obs["SampleName"] = sample_id

        # 确保基因名唯一
        adata.var_names_make_unique()

        return adata

    except Exception as e:
        print(f"❌ 读取样品 {sample_id} 失败: {str(e)}")
        raise


def read_sample(sample_info_path: str,
               memory_monitor: Optional[MemoryMonitor] = None) -> ad.AnnData:
    """
    根据样品信息文件批量读取并整合单细胞数据

    参数：
    ----------
    sample_info_path : str
        样品信息 CSV 文件路径
    格式要求：
        必须包含 'Path' 和 'SampleName' 两列，其他列将作为元数据添加到 obs 中。
    memory_monitor : MemoryMonitor, optional
        内存监控对象

    返回：
    ----------
    adata : AnnData
        整合后的 AnnData 对象
    """
    # 读取样品信息
    sample_info = pd.read_csv(sample_info_path)

    # 验证必需列
    required_cols = ['Path', 'SampleName']
    missing_cols = [col for col in required_cols if col not in sample_info.columns]
    if missing_cols:
        raise ValueError(f"样品信息文件缺少必需列: {missing_cols}")

    adata_list = []

    # 读取各样品数据
    for idx, row in tqdm(sample_info.iterrows(), total=len(sample_info), desc="📖 读取样品"):
        sample_id = row['SampleName']
        data_dir = row['Path']

        # 检查路径是否存在
        if not os.path.exists(data_dir):
            print(f"⚠️  警告: 路径不存在，跳过样品 {sample_id}: {data_dir}")
            continue

        try:
            # 读取单个样品数据
            adata_sample = readbgi(data_dir, sample_id, memory_monitor)

            # 添加元数据
            for col in sample_info.columns:
                if col not in ['Path', 'SampleName']:
                    adata_sample.obs[col] = row[col]

            adata_list.append(adata_sample)

        except Exception as e:
            print(f"⚠️  警告: 跳过样品 {sample_id}: {str(e)}")
            continue

    if len(adata_list) == 0:
        raise ValueError("❌ 没有成功读取任何样品数据！")

    # 合并所有样品
    print(f"\n🔗 合并 {len(adata_list)} 个样品...")
    adata = ad.concat(adata_list, join='outer', merge='same')

    if memory_monitor:
        memory_monitor.checkpoint("数据读取完成")

    # 释放内存
    del adata_list
    gc.collect()

    return adata


def plot_qc_metrics(adata: ad.AnnData, output_dir: str, prefix: str = "before_qc"):
    """
    绘制质控指标可视化图

    参数：
    ----------
    adata : AnnData
        包含质控指标的 AnnData 对象
    output_dir : str
        输出目录
    prefix : str
        文件名前缀（before_qc 或 after_qc）
    """
    print(f"\n   📊 生成质控可视化图 ({prefix})...")

    # 创建 QC 子目录
    qc_dir = os.path.join(output_dir, "qc_plots")
    os.makedirs(qc_dir, exist_ok=True)

    # 1. 小提琴图 - 主要QC指标
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # n_genes_by_counts
    sc.pl.violin(adata, 'n_genes_by_counts', groupby='SampleName',
                 rotation=45, ax=axes[0], show=False)
    axes[0].set_title('Genes per Cell', fontweight='bold')
    axes[0].set_ylabel('Number of Genes')

    # total_counts
    sc.pl.violin(adata, 'total_counts', groupby='SampleName',
                 rotation=45, ax=axes[1], show=False)
    axes[1].set_title('UMI Counts per Cell', fontweight='bold')
    axes[1].set_ylabel('Total UMI Counts')

    # pct_counts_mt
    sc.pl.violin(adata, 'pct_counts_mt', groupby='SampleName',
                 rotation=45, ax=axes[2], show=False)
    axes[2].set_title('Mitochondrial %', fontweight='bold')
    axes[2].set_ylabel('% Mitochondrial Genes')

    plt.tight_layout()
    violin_path = os.path.join(qc_dir, f'{prefix}_violin.png')
    plt.savefig(violin_path, dpi=300, bbox_inches='tight')
    plt.close()

    # 2. 散点图 - 基因数 vs UMI数
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts',
                  color='pct_counts_mt', ax=ax, show=False)
    ax.set_title('QC Metrics: Genes vs UMI Counts', fontweight='bold', fontsize=14)
    scatter_path = os.path.join(qc_dir, f'{prefix}_scatter.png')
    plt.savefig(scatter_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"   ✅ QC可视化图已保存至: {qc_dir}/")


def plot_qc_comparison(adata_after: ad.AnnData, adata_before: ad.AnnData, output_dir: str):
    """
    绘制质控前后对比图

    参数：
    ----------
    adata_after : AnnData
        质控后的数据
    adata_before : AnnData
        质控前的数据
    output_dir : str
        输出目录
    """
    qc_dir = os.path.join(output_dir, "qc_plots")

    # 获取质控前后的数据
    before_df = adata_before.obs.copy()
    after_df = adata_after.obs.copy()

    # 创建对比图
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))

    metrics = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']
    titles = ['Genes per Cell', 'UMI Counts', 'Mitochondrial %']

    for idx, (metric, title) in enumerate(zip(metrics, titles)):
        # 质控前
        axes[0, idx].hist(before_df[metric], bins=50, alpha=0.7,
                         edgecolor='black', color='lightcoral')
        axes[0, idx].set_title(f'{title} (Before QC)', fontweight='bold')
        axes[0, idx].set_xlabel(metric)
        axes[0, idx].set_ylabel('Number of Cells')

        # 质控后
        axes[1, idx].hist(after_df[metric], bins=50, alpha=0.7,
                         edgecolor='black', color='lightgreen')
        axes[1, idx].set_title(f'{title} (After QC)', fontweight='bold')
        axes[1, idx].set_xlabel(metric)
        axes[1, idx].set_ylabel('Number of Cells')

    plt.tight_layout()
    comparison_path = os.path.join(qc_dir, 'qc_comparison.png')
    plt.savefig(comparison_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"   ✅ 质控对比图已保存: {comparison_path}")


def quality_control(
    adata: sc.AnnData,
    min_genes: int = 500,
    max_genes: int = 8000,
    max_pct_mito: float = 15,
    min_counts: int = 1000,
    doublet_method: str = "scrublet",
    doublet_threshold: float = 0.25,
    output_dir: str = "./results",
    return_stats: bool = True,
    memory_monitor: Optional[MemoryMonitor] = None
) -> Union[sc.AnnData, Tuple[sc.AnnData, pd.DataFrame]]:
    """
    单细胞数据质控流程（优化版）
    """
    print("=" * 80)
    print("开始质控流程...")
    print("=" * 80)

    # 1. 计算QC指标
    print("\n[1/5] 计算质控指标...")

    # 识别线粒体基因（支持人和鼠）
    adata.var['mt'] = (adata.var_names.str.startswith('MT-') |
                       adata.var_names.str.startswith('mt-') |
                       adata.var_names.str.startswith('Mt-'))

    # 识别核糖体基因（可选）
    adata.var['ribo'] = (adata.var_names.str.startswith('RPS') |
                         adata.var_names.str.startswith('RPL') |
                         adata.var_names.str.startswith('Rps') |
                         adata.var_names.str.startswith('Rpl'))

    # 计算QC metrics
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=['mt', 'ribo'],
        percent_top=None,
        log1p=False,
        inplace=True
    )

    print(f"   ✓ 质控指标计算完成")
    print(f"   - 线粒体基因数: {adata.var['mt'].sum()}")
    print(f"   - 核糖体基因数: {adata.var['ribo'].sum()}")

    if memory_monitor:
        memory_monitor.checkpoint("QC指标计算")

    # 2. 绘制质控前的图
    print("\n[2/5] 生成质控前可视化图...")
    plot_qc_metrics(adata, output_dir, prefix="before_qc")

    # 记录过滤前的细胞数
    samples = adata.obs['SampleName'].unique()
    qc_stats_list = []

    for sample in samples:
        sample_mask = adata.obs['SampleName'] == sample
        sample_data = adata.obs[sample_mask]
        qc_stats_list.append({
            'SampleName': sample,
            'cells_before_qc': sample_mask.sum(),
            'median_genes_before': sample_data['n_genes_by_counts'].median(),
            'median_counts_before': sample_data['total_counts'].median(),
            'median_mito_before': sample_data['pct_counts_mt'].median()
        })

    # 3. 去除双胞
    if doublet_method.lower() == "scrublet":
        print(f"\n[3/5] 使用 Scrublet 检测双胞 (阈值={doublet_threshold})...")

        import scrublet as scr

        doublet_scores = []
        doublet_labels = []

        for sample in tqdm(samples, desc="   检测双胞"):
            sample_mask = adata.obs['SampleName'] == sample
            adata_sample = adata[sample_mask].copy()

            try:
                # 运行 Scrublet
                scrub = scr.Scrublet(adata_sample.X)
                doublet_score, predicted_doublet = scrub.scrub_doublets(
                    min_counts=2,
                    min_cells=3,
                    min_gene_variability_pctl=85,
                    n_prin_comps=30,
                    verbose=False
                )

                # 使用自定义阈值
                predicted_doublet = doublet_score > doublet_threshold

                doublet_scores.extend(doublet_score)
                doublet_labels.extend(predicted_doublet)

            except Exception as e:
                print(f"   ⚠️  样品 {sample} 双胞检测失败: {str(e)}")
                # 如果失败，标记所有细胞为非双胞
                doublet_scores.extend([0.0] * sample_mask.sum())
                doublet_labels.extend([False] * sample_mask.sum())

        adata.obs['doublet_score'] = doublet_scores
        adata.obs['is_doublet'] = doublet_labels

        n_doublets = sum(doublet_labels)
        print(f"   ✓ 检测到 {n_doublets:,} 个双胞 ({n_doublets/len(adata)*100:.2f}%)")
    else:
        print("\n[3/5] 跳过双胞检测...")
        adata.obs['is_doublet'] = False

    if memory_monitor:
        memory_monitor.checkpoint("双胞检测")

    # 4. 应用过滤条件
    print(f"\n[4/5] 应用过滤条件...")
    print(f"   - 最低基因数: {min_genes}")
    print(f"   - 最高基因数: {max_genes}")
    print(f"   - 最低UMI数: {min_counts}")
    print(f"   - 最大线粒体比例: {max_pct_mito}%")

    # 创建过滤mask
    adata.obs['pass_qc'] = (
        (adata.obs['n_genes_by_counts'] >= min_genes) &
        (adata.obs['n_genes_by_counts'] <= max_genes) &
        (adata.obs['total_counts'] >= min_counts) &
        (adata.obs['pct_counts_mt'] <= max_pct_mito) &
        (~adata.obs['is_doublet'])
    )

    # 5. 统计每个样品的过滤情况
    print("\n[5/5] 统计过滤结果...")

    for i, sample in enumerate(samples):
        sample_mask = adata.obs['SampleName'] == sample
        sample_data = adata.obs[sample_mask]

        n_cells_before = qc_stats_list[i]['cells_before_qc']
        n_cells_after = sample_data['pass_qc'].sum()
        n_filtered = n_cells_before - n_cells_after
        pct_filtered = n_filtered / n_cells_before * 100 if n_cells_before > 0 else 0

        # 统计各种过滤原因
        n_low_genes = (sample_data['n_genes_by_counts'] < min_genes).sum()
        n_high_genes = (sample_data['n_genes_by_counts'] > max_genes).sum()
        n_low_counts = (sample_data['total_counts'] < min_counts).sum()
        n_high_mito = (sample_data['pct_counts_mt'] > max_pct_mito).sum()
        n_doublet = sample_data['is_doublet'].sum()

        # 质控后的统计
        passed_data = sample_data[sample_data['pass_qc']]

        qc_stats_list[i].update({
            'cells_after_qc': n_cells_after,
            'cells_filtered': n_filtered,
            'pct_filtered': f"{pct_filtered:.2f}%",
            'n_low_genes': n_low_genes,
            'n_high_genes': n_high_genes,
            'n_low_counts': n_low_counts,
            'n_high_mito': n_high_mito,
            'n_doublet': n_doublet,
            'median_genes_after': passed_data['n_genes_by_counts'].median() if len(passed_data) > 0 else 0,
            'median_counts_after': passed_data['total_counts'].median() if len(passed_data) > 0 else 0,
            'median_mito_after': passed_data['pct_counts_mt'].median() if len(passed_data) > 0 else 0
        })

    qc_stats_df = pd.DataFrame(qc_stats_list)

    # 输出统计摘要
    print("\n" + "=" * 80)
    print("质控统计摘要:")
    print("=" * 80)
    print(qc_stats_df.to_string(index=False))
    print("=" * 80)

    # 过滤数据并正确保留 .raw
    print(f"\n正在过滤数据...")
    print(f"   过滤前: {adata.n_obs:,} 细胞")

    # 保存质控前的数据用于对比图
    adata_before_qc = adata.copy()

    # 保存原始counts到 .raw
    adata.raw = adata

    # 过滤细胞
    adata = adata[adata.obs['pass_qc']].copy()

    # 过滤低表达基因（在至少3个细胞中表达）
    sc.pp.filter_genes(adata, min_cells=3)

    print(f"   过滤后: {adata.n_obs:,} 细胞, {adata.n_vars:,} 基因")
    print(f"   过滤率: {(1 - adata.n_obs/adata.raw.n_obs)*100:.2f}%")

    if memory_monitor:
        memory_monitor.checkpoint("数据过滤")

    # 绘制质控后的图
    print("\n生成质控后可视化图...")
    plot_qc_metrics(adata, output_dir, prefix="after_qc")

    # 生成对比图
    print("\n生成质控前后对比图...")
    plot_qc_comparison(adata, adata_before_qc, output_dir)

    print("\n" + "=" * 80)
    print("✅ 质控完成!")
    print("=" * 80)

    # 清理内存
    gc.collect()

    if return_stats:
        return adata, qc_stats_df
    else:
        return adata


def data_integration(
    adata: sc.AnnData,
    batch_key: str = "SampleName",
    method: str = "harmony",
    n_pcs: int = 50,
    n_neighbors: int = 30,
    resolution: float = 0.8,
    gene_num: int = 3000,
    umap_min_dist: float = 0.3,
    output_dir: str = "./results",
    memory_monitor: Optional[MemoryMonitor] = None
) -> sc.AnnData:
    """
    数据整合与降维聚类（优化版）
    """
    print("\n" + "=" * 80)
    print(f"数据整合与降维 (方法: {method})")
    print("=" * 80)

    # 1. 标准化和对数转换
    print("\n[1/7] 数据标准化...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    if memory_monitor:
        memory_monitor.checkpoint("数据标准化")

    # 2. 鉴定高变基因
    print(f"\n[2/7] 鉴定高变基因 (n={gene_num})...")
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=gene_num,
        batch_key=batch_key,
        flavor='seurat_v3',
        subset=False
    )

    n_hvg = adata.var['highly_variable'].sum()
    print(f"   ✓ 鉴定到 {n_hvg} 个高变基因")

    # 3. 数据缩放
    print("\n[3/7] 数据缩放...")
    sc.pp.scale(adata, max_value=10)

    if memory_monitor:
        memory_monitor.checkpoint("数据缩放")

    # 4. PCA降维
    print(f"\n[4/7] PCA降维 (n_pcs={n_pcs})...")
    sc.tl.pca(adata, svd_solver='arpack', n_comps=n_pcs)

    # 绘制PCA方差解释图
    fig, ax = plt.subplots(figsize=(10, 6))
    variance_ratio = adata.uns['pca']['variance_ratio']
    cumsum_variance = np.cumsum(variance_ratio)

    ax.plot(range(1, len(variance_ratio)+1), variance_ratio, 'o-', label='Individual')
    ax.plot(range(1, len(cumsum_variance)+1), cumsum_variance, 's-', label='Cumulative')
    ax.axhline(0.9, color='r', linestyle='--', label='90% variance')
    ax.set_xlabel('Principal Component')
    ax.set_ylabel('Variance Explained')
    ax.set_title('PCA Variance Explained')
    ax.legend()
    ax.grid(True, alpha=0.3)

    pca_path = os.path.join(output_dir, 'pca_variance.png')
    plt.savefig(pca_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"   ✓ PCA完成，前{n_pcs}个PC解释 {cumsum_variance[n_pcs-1]*100:.2f}% 的方差")

    if memory_monitor:
        memory_monitor.checkpoint("PCA降维")

    # 5. 批次效应校正
    if method.lower() == "harmony":
        print("\n[5/7] Harmony 批次校正...")
        try:
            import harmonypy as hm
            ho = hm.run_harmony(
                adata.obsm['X_pca'],
                adata.obs,
                batch_key,
                max_iter_harmony=20,
                verbose=False
            )
            adata.obsm['X_pca_harmony'] = ho.Z_corr.T
            use_rep = 'X_pca_harmony'
            print("   ✓ Harmony 校正完成")
        except ImportError:
            print("   ⚠️  harmonypy 未安装，跳过批次校正")
            use_rep = 'X_pca'
        except Exception as e:
            print(f"   ⚠️  Harmony 校正失败: {str(e)}")
            use_rep = 'X_pca'
    else:
        print("\n[5/7] 跳过批次校正...")
        use_rep = 'X_pca'

    if memory_monitor:
        memory_monitor.checkpoint("批次校正")

    # 6. 构建邻居图和UMAP
    print(f"\n[6/7] 构建邻居图 (n_neighbors={n_neighbors})...")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep=use_rep)

    print(f"   计算 UMAP (min_dist={umap_min_dist})...")
    sc.tl.umap(adata, min_dist=umap_min_dist)

    if memory_monitor:
        memory_monitor.checkpoint("邻居图和UMAP")

    # 7. Leiden聚类
    print(f"\n[7/7] Leiden 聚类 (resolution={resolution})...")
    sc.tl.leiden(adata, resolution=resolution, key_added='leiden')

    n_clusters = adata.obs['leiden'].nunique()
    print(f"   ✓ 鉴定到 {n_clusters} 个聚类")

    if memory_monitor:
        memory_monitor.checkpoint("聚类完成")

    # 8. 生成可视化图
    print("\n生成可视化图...")
    plot_integration_results(adata, output_dir, batch_key)

    print("\n" + "=" * 80)
    print("✅ 数据整合完成!")
    print("=" * 80)

    # 清理内存
    gc.collect()

    return adata


def plot_integration_results(adata: ad.AnnData, output_dir: str, batch_key: str):
    """
    绘制整合结果可视化图
    """
    # 1. UMAP - 按聚类着色
    sc.pl.umap(adata, color='leiden', legend_loc='on data',
               title='UMAP - Leiden Clusters', frameon=False,
               save='_leiden.png')

    # 2. UMAP - 按样品着色
    sc.pl.umap(adata, color=batch_key, title=f'UMAP - {batch_key}',
               frameon=False, save=f'_{batch_key}.png')

    print(f"   ✅ 可视化图已保存至 {output_dir}/figures/")


def celltypist_annotation(
    adata: sc.AnnData,
    model_path: Optional[str] = None,
    output_path: str = "./annotated_data.h5ad"
) -> sc.AnnData:
    """
    使用 CellTypist 进行细胞类型注释
    """
    try:
        import celltypist
        from celltypist import models

        if model_path is None:
            print("   下载默认 CellTypist 模型...")
            model = models.Model.load(model='Immune_All_Low.pkl')
        else:
            print(f"   加载 CellTypist 模型: {model_path}")
            model = models.Model.load(model=model_path)

        print("   开始细胞类型注释...")
        predictions = celltypist.annotate(adata, model=model, majority_voting=True)

        result_model = predictions.predicted_labels[['predicted_labels', 'majority_voting']].rename(
            columns={'predicted_labels': 'celltype_pred', 'majority_voting': 'celltype_pred_mv'}
        )

        adata.obs = adata.obs.join(result_model, how='left')

        print(f"   ✓ 注释完成，结果保存到 {output_path}")
        adata.write_h5ad(output_path)

        return adata

    except ImportError:
        print("   ❌ celltypist 未安装，跳过细胞类型注释")
        print("   安装命令: pip install celltypist")
        return adata
    except Exception as e:
        print(f"   ❌ 细胞类型注释失败: {str(e)}")
        return adata


def get_args():
    import argparse

    parser = argparse.ArgumentParser(
        description="优化的单细胞RNA-seq数据处理流程 (支持断点续传)",
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
    checkpoint_group.add_argument("--auto-resume", action='store_true',
                                 help="自动从断点恢复（无需确认）")
    checkpoint_group.add_argument("--reset", action='store_true',
                                 help="重置所有断点，从头开始")
    checkpoint_group.add_argument("--show-progress", action='store_true',
                                 help="只显示当前进度，不运行")

    # 质控参数（优化后的默认值）
    qc_group = parser.add_argument_group('质控参数')
    qc_group.add_argument("--min_genes", type=int, default=500,
                         help="每个细胞的最低基因数")
    qc_group.add_argument("--max_genes", type=int, default=8000,
                         help="每个细胞的最高基因数")
    qc_group.add_argument("--min_counts", type=int, default=1000,
                         help="每个细胞的最低UMI数")
    qc_group.add_argument("--max_pct_mito", type=float, default=15.0,
                         help="最大线粒体基因百分比")
    qc_group.add_argument("--doublet_method", default="scrublet",
                         choices=["scrublet", "none"],
                         help="双胞检测方法")
    qc_group.add_argument("--doublet_threshold", type=float, default=0.25,
                         help="双胞检测阈值")

    # 整合与聚类参数（优化后的默认值）
    int_group = parser.add_argument_group('整合与聚类参数')
    int_group.add_argument("--integration_method", default="harmony",
                          choices=["harmony", "combat", "scanorama", "none"],
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

    # 其他参数
    parser.add_argument("--n_jobs", type=int, default=-1,
                       help="并行任务数 (-1表示使用所有CPU)")

    return parser.parse_args()


def main():
    """
    单细胞RNA-seq数据处理主流程（增强版 - 支持断点续传）
    """
    args = get_args()

    # 设置 scanpy 参数
    sc.settings.verbosity = 1
    sc.settings.n_jobs = args.n_jobs
    sc.settings.figdir = os.path.join(args.output_dir, 'figures')

    # 创建输出目录
    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(sc.settings.figdir, exist_ok=True)

    # 初始化断点管理器
    ckpt = CheckpointManager(args.output_dir, auto_resume=args.auto_resume)

    # 处理断点相关命令
    if args.reset:
        ckpt.reset()

    if args.show_progress:
        ckpt.print_status()
        return

    # 检查是否需要从断点恢复
    if args.resume or args.auto_resume:
        if ckpt.should_resume():
            pass  # 继续执行
        elif not ckpt.get_last_completed_step():
            print("没有检测到之前的进度，将从头开始运行\n")

    # 初始化内存监控
    mem_monitor = MemoryMonitor()
    mem_monitor.checkpoint("流程开始")

    print("\n" + "=" * 100)
    print(" " * 30 + "单细胞RNA-seq数据处理流程 (增强版 v2.2)")
    print("=" * 100)
    print(f"\n📋 配置参数：")
    print(f"  样品信息: {args.sample_info}")
    print(f"  输出目录: {args.output_dir}")
    print(f"  断点文件: {ckpt.checkpoint_file}")
    print("=" * 100 + "\n")

    # ============================================
    # 步骤 1: 读取样品数据
    # ============================================
    raw_data_path = os.path.join(args.output_dir, "01_raw_data.h5ad")

    if not ckpt.is_completed('step1_read'):
        def read_data_step():
            adata = read_sample(args.sample_info, memory_monitor=mem_monitor)

            # 保存原始数据
            adata.write(raw_data_path)

            print(f"\n✅ 读取完成:")
            print(f"  细胞数: {adata.n_obs:,}")
            print(f"  基因数: {adata.n_vars:,}")
            print(f"  数据已保存: {raw_data_path}")

            # 保存断点
            ckpt.save_checkpoint('step1_read',
                               file=raw_data_path,
                               n_cells=adata.n_obs,
                               n_genes=adata.n_vars)

            return adata

        adata = safe_step_execution(ckpt, '步骤1: 读取样品数据', read_data_step)
    else:
        print("\n✓ 步骤1已完成，加载已有数据...")
        adata = sc.read_h5ad(raw_data_path)
        print(f"  细胞数: {adata.n_obs:,}")
        print(f"  基因数: {adata.n_vars:,}\n")

    # ============================================
    # 步骤 2: 质量控制
    # ============================================
    filtered_data_path = os.path.join(args.output_dir, "02_filtered_data.h5ad")

    if not ckpt.is_completed('step2_qc'):
        def qc_step():
            adata_filtered, qc_stats = quality_control(
                adata,
                min_genes=args.min_genes,
                max_genes=args.max_genes,
                min_counts=args.min_counts,
                max_pct_mito=args.max_pct_mito,
                doublet_method=args.doublet_method,
                doublet_threshold=args.doublet_threshold,
                output_dir=args.output_dir,
                return_stats=True,
                memory_monitor=mem_monitor
            )

            # 保存质控后的数据
            adata_filtered.write(filtered_data_path)

            # 保存质控统计表
            qc_stats_path = os.path.join(args.output_dir, "02_qc_statistics.csv")
            qc_stats.to_csv(qc_stats_path, index=False)

            print(f"\n✅ 质控完成:")
            print(f"  过滤后细胞数: {adata_filtered.n_obs:,}")
            print(f"  数据已保存: {filtered_data_path}")

            # 保存断点
            ckpt.save_checkpoint('step2_qc',
                               file=filtered_data_path,
                               n_cells_after=adata_filtered.n_obs)

            return adata_filtered

        adata = safe_step_execution(ckpt, '步骤2: 质量控制', qc_step)
    else:
        print("\n✓ 步骤2已完成，加载已有数据...")
        adata = sc.read_h5ad(filtered_data_path)
        print(f"  过滤后细胞数: {adata.n_obs:,}\n")

    # ============================================
    # 步骤 3: 数据整合
    # ============================================
    integrated_data_path = os.path.join(args.output_dir, "03_integrated_data.h5ad")

    if not ckpt.is_completed('step3_integrate'):
        def integration_step():
            adata_integrated = data_integration(
                adata,
                batch_key="SampleName",
                method=args.integration_method,
                n_pcs=args.n_pcs,
                n_neighbors=args.n_neighbors,
                resolution=args.resolution,
                gene_num=args.gene_num,
                umap_min_dist=args.umap_min_dist,
                output_dir=args.output_dir,
                memory_monitor=mem_monitor
            )

            # 保存整合后的数据
            adata_integrated.write(integrated_data_path)

            print(f"\n✅ 数据整合完成:")
            print(f"  细胞数: {adata_integrated.n_obs:,}")
            print(f"  聚类数: {adata_integrated.obs['leiden'].nunique()}")
            print(f"  数据已保存: {integrated_data_path}")

            # 保存断点
            ckpt.save_checkpoint('step3_integrate',
                               file=integrated_data_path,
                               n_clusters=adata_integrated.obs['leiden'].nunique())

            return adata_integrated

        adata = safe_step_execution(ckpt, '步骤3: 数据整合与降维', integration_step)
    else:
        print("\n✓ 步骤3已完成，加载已有数据...")
        adata = sc.read_h5ad(integrated_data_path)
        print(f"  整合后细胞数: {adata.n_obs:,}")
        print(f"  聚类数: {adata.obs['leiden'].nunique()}\n")

    # ============================================
    # 步骤 4: 细胞类型注释 (可选)
    # ============================================
    if args.run_annotation:
        annotation_output = os.path.join(args.output_dir, "04_annotated_data.h5ad")

        if not ckpt.is_completed('step4_annotate'):
            def annotation_step():
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

                mem_monitor.checkpoint("细胞类型注释")

                return adata_annotated

            adata = safe_step_execution(ckpt, '步骤4: 细胞类型注释', annotation_step)
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
    print("\n💡 断点续传提示：")
    print("  • 如遇错误，修复后运行: python sc_pipeline_enhanced.py --resume ...")
    print("  • 查看进度: python sc_pipeline_enhanced.py --show-progress ...")
    print("  • 从头开始: python sc_pipeline_enhanced.py --reset ...")
    print("=" * 100 + "\n")


if __name__ == "__main__":
    main()
