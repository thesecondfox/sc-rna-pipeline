# 断点续传功能使用指南

## 概述

增强版单细胞RNA-seq分析流程支持完整的断点续传功能。当流程因错误中断时，可以修复问题后从中断点继续运行，无需从头开始。

## 核心功能

### 1. 自动断点保存

流程会在每个步骤完成后自动保存断点：

- **步骤1**: 数据读取完成后保存 `01_raw_data.h5ad`
- **步骤2**: 质量控制完成后保存 `02_filtered_data.h5ad`
- **步骤3**: 数据整合完成后保存 `03_integrated_data.h5ad`
- **步骤4**: 细胞注释完成后保存 `04_annotated_data.h5ad`

### 2. 数据完整性检查

- 自动验证断点文件是否存在
- 检查文件大小是否匹配
- 标记损坏的断点为无效

### 3. 错误处理机制

- 捕获所有异常并记录详细信息
- 保存错误日志到 `.error_log.txt`
- 在断点文件中标记失败的步骤

### 4. 交互式恢复

- 检测到断点时提示用户确认
- 显示详细的进度信息
- 支持自动恢复模式

## 使用方法

### 基础运行（首次运行）

```bash
python sc_pipeline_enhanced.py \
    --sample_info samples.csv \
    --output_dir ./results \
    --min_genes 500 \
    --max_genes 8000 \
    --integration_method harmony
```

### 从断点恢复（交互式）

如果流程中断，使用 `--resume` 参数：

```bash
python sc_pipeline_enhanced.py \
    --sample_info samples.csv \
    --output_dir ./results \
    --resume
```

系统会显示进度并询问是否继续：

```
检测到之前的进度:
================================================================================
当前分析进度
================================================================================
✅ 步骤1: 数据读取
   完成时间: 2025-10-28 10:30:15
   文件: 01_raw_data.h5ad (1250.5 MB)
   细胞数: 45,230
   基因数: 33,538

✅ 步骤2: 质量控制
   完成时间: 2025-10-28 10:45:22
   文件: 02_filtered_data.h5ad (980.3 MB)
   细胞数: 38,567

❌ 步骤3: 数据整合 - 失败
   错误: Memory allocation error
   时间: 2025-10-28 10:50:18

⏳ 步骤4: 细胞注释 - 待执行
================================================================================

是否从上次断点继续运行？(y/n):
```

### 自动恢复（无需确认）

适用于自动化脚本：

```bash
python sc_pipeline_enhanced.py \
    --sample_info samples.csv \
    --output_dir ./results \
    --auto-resume
```

### 查看当前进度

不运行流程，只查看进度：

```bash
python sc_pipeline_enhanced.py \
    --sample_info samples.csv \
    --output_dir ./results \
    --show-progress
```

### 重置断点（从头开始）

删除所有断点，重新运行：

```bash
python sc_pipeline_enhanced.py \
    --sample_info samples.csv \
    --output_dir ./results \
    --reset
```

## 断点文件说明

### `.checkpoint.json`

存储所有断点信息的JSON文件：

```json
{
  "step1_read": {
    "completed": true,
    "timestamp": "2025-10-28T10:30:15.123456",
    "valid": true,
    "file": "./results/01_raw_data.h5ad",
    "file_size": 1311551488,
    "n_cells": 45230,
    "n_genes": 33538
  },
  "step2_qc": {
    "completed": true,
    "timestamp": "2025-10-28T10:45:22.789012",
    "valid": true,
    "file": "./results/02_filtered_data.h5ad",
    "file_size": 1027995648,
    "n_cells_after": 38567
  },
  "step3_integrate": {
    "failed": true,
    "timestamp": "2025-10-28T10:50:18.456789",
    "error_type": "MemoryError",
    "error_message": "Memory allocation error",
    "traceback": "..."
  }
}
```

### `.error_log.txt`

详细的错误日志文件，包含：

- 错误步骤和时间
- 错误类型和消息
- 完整的堆栈追踪

## 常见场景

### 场景1: 内存不足导致中断

**问题**: 在数据整合步骤因内存不足失败

**解决方案**:
1. 释放内存或增加系统内存
2. 降低参数（如 `--gene_num 2000` 或 `--n_pcs 30`）
3. 使用 `--resume` 继续运行

```bash
# 修改参数后继续
python sc_pipeline_enhanced.py \
    --sample_info samples.csv \
    --output_dir ./results \
    --gene_num 2000 \
    --n_pcs 30 \
    --resume
```

### 场景2: 网络问题导致数据读取失败

**问题**: 某个样品文件读取失败

**解决方案**:
1. 检查并修复文件路径或下载缺失文件
2. 更新 `samples.csv`
3. 使用 `--reset` 重新开始（因为步骤1未完成）

```bash
python sc_pipeline_enhanced.py \
    --sample_info samples.csv \
    --output_dir ./results \
    --reset
```

### 场景3: 质控参数不合适

**问题**: 完成质控后发现参数过于严格，过滤掉太多细胞

**解决方案**:
1. 调整质控参数
2. 删除步骤2的断点，重新运行质控

```bash
# 删除 02_filtered_data.h5ad 和修改 .checkpoint.json
rm ./results/02_filtered_data.h5ad

# 重新运行，会自动跳过步骤1
python sc_pipeline_enhanced.py \
    --sample_info samples.csv \
    --output_dir ./results \
    --min_genes 300 \
    --max_pct_mito 20 \
    --resume
```

### 场景4: 测试不同的整合方法

**问题**: 想测试不同的批次校正方法

**解决方案**:
1. 创建不同的输出目录
2. 从步骤2的结果开始

```bash
# 方法1: Harmony
python sc_pipeline_enhanced.py \
    --sample_info samples.csv \
    --output_dir ./results_harmony \
    --integration_method harmony \
    --resume

# 方法2: Combat
python sc_pipeline_enhanced.py \
    --sample_info samples.csv \
    --output_dir ./results_combat \
    --integration_method combat \
    --resume
```

## 进度恢复逻辑

流程会按以下逻辑判断从哪一步继续：

```
检查 step4_annotate → 已完成？跳过
检查 step3_integrate → 已完成？跳过，否则执行
检查 step2_qc → 已完成？跳过，否则执行
检查 step1_read → 已完成？跳过，否则执行
```

## 最佳实践

### 1. 定期备份

虽然有断点续传，仍建议定期备份重要数据文件：

```bash
# 备份中间结果
cp -r ./results ./results_backup_$(date +%Y%m%d)
```

### 2. 使用独立的输出目录

不同的分析使用不同的输出目录：

```bash
python sc_pipeline_enhanced.py \
    --sample_info samples.csv \
    --output_dir ./results_v1 \
    --resolution 0.8

python sc_pipeline_enhanced.py \
    --sample_info samples.csv \
    --output_dir ./results_v2 \
    --resolution 1.2
```

### 3. 查看错误日志

遇到错误时，查看详细日志：

```bash
cat ./results/.error_log.txt
```

### 4. 验证数据完整性

使用Python验证保存的数据：

```python
import scanpy as sc

# 读取数据
adata = sc.read_h5ad('./results/02_filtered_data.h5ad')

# 检查数据
print(f"细胞数: {adata.n_obs:,}")
print(f"基因数: {adata.n_vars:,}")
print(f"包含 .raw: {adata.raw is not None}")
```

## 故障排查

### 断点文件损坏

如果 `.checkpoint.json` 损坏：

```bash
# 删除断点文件
rm ./results/.checkpoint.json

# 但保留数据文件，系统会自动重建断点
python sc_pipeline_enhanced.py \
    --sample_info samples.csv \
    --output_dir ./results \
    --resume
```

### 数据文件损坏

如果某个步骤的数据文件损坏：

1. 删除该步骤及后续步骤的数据文件
2. 编辑 `.checkpoint.json`，移除对应步骤
3. 使用 `--resume` 继续

```bash
# 例如：步骤3的数据损坏
rm ./results/03_integrated_data.h5ad
rm ./results/04_annotated_data.h5ad

# 手动编辑 .checkpoint.json 或直接删除
rm ./results/.checkpoint.json

python sc_pipeline_enhanced.py \
    --sample_info samples.csv \
    --output_dir ./results \
    --resume
```

## 性能优化建议

### 内存管理

1. **监控内存使用**: 查看 `memory_usage.csv` 和 `memory_usage.png`
2. **调整参数**: 降低 `--gene_num` 和 `--n_pcs` 可减少内存占用
3. **分批处理**: 如果样品太多，考虑分批处理

### 并行处理

使用 `--n_jobs` 参数控制并行度：

```bash
# 使用所有CPU核心
python sc_pipeline_enhanced.py --n_jobs -1 ...

# 使用4个核心
python sc_pipeline_enhanced.py --n_jobs 4 ...
```

## 总结

增强版流程的断点续传功能提供了：

- ✅ **自动保存**: 每步完成自动保存
- ✅ **完整性验证**: 检查数据文件完整性
- ✅ **错误处理**: 详细的错误日志
- ✅ **灵活恢复**: 交互式或自动恢复
- ✅ **进度查看**: 随时查看当前状态

使您的分析流程更加可靠和高效！
