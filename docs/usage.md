# 使用指南

本文档详细介绍如何使用单细胞 RNA-seq 分析流程。

## 目录

- [基础用法](#基础用法)
- [输入文件准备](#输入文件准备)
- [参数详解](#参数详解)
- [高级用法](#高级用法)
- [HPC 使用](#hpc-使用)
- [结果解读](#结果解读)

## 基础用法

### 最简单的运行方式
```bash
python sc_pipeline.py --sample_info samples.csv --output_dir ./results
```

这将使用所有默认参数运行完整的分析流程。

### 查看帮助
```bash
python sc_pipeline.py --help
```

## 输入文件准备

### 样品信息文件 (samples.csv)

这是一个 CSV 格式的文件，包含所有样品的信息。

#### 必需列

| 列名 | 说明 | 示例 |
|------|------|------|
| Path | 数据目录路径 | /data/sample1 |
| SampleName | 样品唯一ID | Sample1 |

#### 可选列

任何其他列都会被保留为元数据，例如：Group, Stage, Region, Batch, Gender, Age 等。

#### 示例
```csv
Path,SampleName,Group,Stage,Region
/data/sample1,Sample1,Control,Adult,Colon
/data/sample2,Sample2,Treatment,Adult,Colon
/data/sample3,Sample3,Control,Fetal,Ileum
```

### 数据目录结构

每个 `Path` 指向的目录应包含以下文件：
```
sample1/
├── matrix.mtx(.gz)      # 表达矩阵（稀疏格式）
├── barcodes.tsv(.gz)    # 细胞条码
└── features.tsv(.gz)    # 基因名称
```

支持压缩 (.gz) 和未压缩格式。

## 参数详解

### 质控参数

#### --min_genes (默认: 200)
每个细胞检测到的最低基因数。
```bash
--min_genes 200   # 宽松
--min_genes 500   # 严格
```

#### --max_genes (默认: 6000)
每个细胞检测到的最高基因数。
```bash
--max_genes 6000   # 标准
--max_genes 8000   # 宽松
```

#### --max_pct_mito (默认: 20)
线粒体基因占比的最大百分比。
```bash
--max_pct_mito 20   # 标准
--max_pct_mito 15   # 严格
```

#### --doublet_method (默认: scrublet)
双胞检测方法。
```bash
--doublet_method scrublet   # 使用 Scrublet
--doublet_method none       # 跳过双胞检测
```

#### --doublet_threshold (默认: 0.25)
Scrublet 双胞评分阈值。

### 整合参数

#### --integration_method (默认: harmony)
批次效应校正方法。
```bash
--integration_method harmony   # Harmony（推荐）
--integration_method combat    # ComBat
--integration_method none      # 不进行批次校正
```

#### --n_pcs (默认: 30)
主成分分析保留的主成分数量。

#### --n_neighbors (默认: 10)
构建 KNN 图时的邻居数量。

#### --resolution (默认: 1.1)
Leiden 聚类的分辨率参数。
```bash
--resolution 0.5   # 更少的 cluster
--resolution 1.1   # 标准
--resolution 2.0   # 更多的 cluster
```

#### --gene_num (默认: 2000)
选择的高变基因数量。

#### --umap_min_dist (默认: 0.5)
UMAP 的最小距离参数。

### 注释参数

#### --run_annotation
是否运行细胞类型注释（CellTypist）。

#### --celltypist_model
CellTypist 模型路径。
```bash
--celltypist_model /path/to/model.pkl
```

### 输出参数

#### --output_dir (默认: ./results)
输出目录路径。

## 高级用法

### 分步运行
```python
import scanpy as sc
from sc_pipeline import read_sample, quality_control, data_integration

# 步骤 1: 读取数据
adata = read_sample("samples.csv")
adata.write("raw_data.h5ad")

# 步骤 2: 质控
adata_filtered, qc_stats = quality_control(adata)
adata_filtered.write("filtered_data.h5ad")

# 步骤 3: 整合
adata_integrated = data_integration(adata_filtered)
adata_integrated.write("integrated_data.h5ad")
```

### 批量处理
```bash
#!/bin/bash

projects=("project1" "project2" "project3")

for proj in "${projects[@]}"; do
    python sc_pipeline.py \
        --sample_info ${proj}/samples.csv \
        --output_dir ${proj}/results
done
```

## HPC 使用

### 使用 LSF

编辑 `submit_job.sh` 然后提交:
```bash
./submit_job.sh
```

### 使用 Slurm

创建 `submit_slurm.sh`:
```bash
#!/bin/bash
#SBATCH --job-name=sc-pipeline
#SBATCH --output=sc_pipeline.out
#SBATCH --error=sc_pipeline.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G

python sc_pipeline.py \
    --sample_info samples.csv \
    --output_dir ./results
```

提交:
```bash
sbatch submit_slurm.sh
```

### 监控任务

**LSF**:
```bash
bjobs              # 查看所有任务
bjobs -l <job_id>  # 查看详细信息
bkill <job_id>     # 取消任务
```

**Slurm**:
```bash
squeue -u $USER    # 查看任务
scancel <job_id>   # 取消任务
```

## 结果解读

### 输出文件

| 文件 | 说明 |
|------|------|
| `raw_data.h5ad` | 原始数据 |
| `filtered_data.h5ad` | 质控后数据 |
| `integrated_data.h5ad` | 整合后数据 |
| `qc_statistics.csv` | 质控统计 |
| `umap_*.png` | UMAP 可视化 |
| `heterogeneity_*.png` | 异质性分析 |

### 读取结果
```python
import scanpy as sc
import pandas as pd

# 读取数据
adata = sc.read_h5ad("results/integrated_data.h5ad")

# 查看基本信息
print(f"细胞数: {adata.n_obs}")
print(f"基因数: {adata.n_vars}")

# 查看质控统计
qc_stats = pd.read_csv("results/qc_statistics.csv")
print(qc_stats)
```

## 最佳实践

1. **先用小数据集测试**: 用 1-2 个样品测试参数
2. **查看质控报告**: 仔细检查过滤统计
3. **检查批次效应**: 看整合前后的 UMAP 对比
4. **验证聚类**: 使用已知 marker 基因验证
5. **保存中间文件**: 便于重新分析

## 下一步

- 查看 [完整教程](tutorial.md) 了解详细示例
- 查看 [FAQ](faq.md) 了解常见问题
