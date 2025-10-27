# 断点续传功能说明

## 概述

sc-rna-pipeline v2.0 支持**流程断点续传**功能，当分析流程因为报错、中断或其他原因停止后，可以从上次中断的地方继续运行，无需重新开始。

## 功能特性

### 1. 自动保存断点
- ✅ 每个主要步骤完成后自动保存进度
- ✅ 保存中间结果文件路径
- ✅ 记录完成时间和关键统计信息
- ✅ 断点文件：`.checkpoint.json`

### 2. 智能恢复
- ✅ 自动检测已完成的步骤
- ✅ 跳过已完成的步骤，直接加载结果
- ✅ 支持从任意断点继续
- ✅ 支持完全重置，从头开始

### 3. 灵活控制
- ✅ 查看当前进度
- ✅ 跳过到指定步骤
- ✅ 手动重置断点
- ✅ 自动内存管理

## 使用方法

### 基本用法

#### 1. 正常运行（第一次）

```bash
sc-pipeline \
    --sample_info samples.csv \
    --output_dir ./results \
    --min_genes 500 \
    --max_genes 8000
```

如果流程中断（报错、手动停止等），会自动保存断点到 `./results/.checkpoint.json`

#### 2. 从断点继续运行

```bash
sc-pipeline \
    --sample_info samples.csv \
    --output_dir ./results \
    --resume
```

使用 `--resume` 参数，流程会：
1. 读取断点文件
2. 显示当前进度
3. 跳过已完成的步骤
4. 从中断点继续运行

#### 3. 查看当前进度

```bash
sc-pipeline \
    --sample_info samples.csv \
    --output_dir ./results \
    --show-progress
```

输出示例：
```
================================================================================
当前分析进度
================================================================================
✅ 步骤1: 数据读取
   完成时间: 2025-01-15 10:30:45
   文件: 01_raw_data.h5ad (1250.3 MB)
✅ 步骤2: 质量控制
   完成时间: 2025-01-15 10:45:20
   文件: 02_filtered_data.h5ad (980.5 MB)
⏳ 步骤3: 数据整合
⏳ 步骤4: 细胞注释
================================================================================
```

#### 4. 重置所有断点，从头开始

```bash
sc-pipeline \
    --sample_info samples.csv \
    --output_dir ./results \
    --reset
```

这会删除断点文件，从第一步重新开始。

#### 5. 跳过到指定步骤

```bash
# 跳过到质控步骤
sc-pipeline \
    --sample_info samples.csv \
    --output_dir ./results \
    --skip-to qc

# 跳过到整合步骤
sc-pipeline \
    --sample_info samples.csv \
    --output_dir ./results \
    --skip-to integrate

# 跳过到注释步骤
sc-pipeline \
    --sample_info samples.csv \
    --output_dir ./results \
    --skip-to annotate
```

## 工作流程示例

### 场景1：流程在质控步骤报错

```bash
# 第一次运行
$ sc-pipeline --sample_info samples.csv --output_dir ./results

# 输出：
🔹🔹🔹... 步骤 1/4: 读取样品数据
✅ 读取完成
   💾 断点已保存: step1_read

🔹🔹🔹... 步骤 2/4: 质量控制
❌ 质控失败: XXX错误

# 修复问题后，从断点继续
$ sc-pipeline --sample_info samples.csv --output_dir ./results --resume

# 输出：
================================================================================
当前分析进度
================================================================================
✅ 步骤1: 数据读取
⏳ 步骤2: 质量控制
⏳ 步骤3: 数据整合
⏳ 步骤4: 细胞注释
================================================================================
🔄 从上次断点继续运行...

✓ 步骤1已完成，加载已有数据...
  细胞数: 50,000

🔹🔹🔹... 步骤 2/4: 质量控制
# 继续运行质控...
```

### 场景2：手动中断后继续

```bash
# 运行中按 Ctrl+C 中断
$ sc-pipeline --sample_info samples.csv --output_dir ./results
^C

# 查看进度
$ sc-pipeline --sample_info samples.csv --output_dir ./results --show-progress
✅ 步骤1: 数据读取
✅ 步骤2: 质量控制
⏳ 步骤3: 数据整合

# 继续运行
$ sc-pipeline --sample_info samples.csv --output_dir ./results --resume
```

### 场景3：修改参数后重新运行

```bash
# 第一次运行，发现参数不合适
$ sc-pipeline \
    --sample_info samples.csv \
    --output_dir ./results \
    --min_genes 200  # 太宽松

# 重置断点，使用新参数重新运行
$ sc-pipeline \
    --sample_info samples.csv \
    --output_dir ./results \
    --min_genes 500 \  # 更严格
    --reset
```

## 断点文件结构

`.checkpoint.json` 示例：

```json
{
  "step1_read": {
    "completed": true,
    "timestamp": "2025-01-15T10:30:45.123456",
    "file": "./results/01_raw_data.h5ad",
    "n_cells": 50000,
    "n_genes": 30000
  },
  "step2_qc": {
    "completed": true,
    "timestamp": "2025-01-15T10:45:20.654321",
    "file": "./results/02_filtered_data.h5ad",
    "n_cells_after": 45000
  },
  "step3_integrate": {
    "completed": false
  }
}
```

## Python API 使用

```python
from sc_pipeline import CheckpointManager, MemoryMonitor
from sc_pipeline.io import read_sample
from sc_pipeline.preprocessing import quality_control

# 初始化断点管理器
ckpt = CheckpointManager('./results')

# 查看进度
ckpt.print_status()

# 步骤1：读取数据
if not ckpt.is_completed('step1_read'):
    adata = read_sample('samples.csv')
    adata.write('./results/01_raw_data.h5ad')
    ckpt.save_checkpoint('step1_read',
                        file='./results/01_raw_data.h5ad',
                        n_cells=adata.n_obs)
else:
    print("步骤1已完成，加载数据...")
    adata = sc.read_h5ad('./results/01_raw_data.h5ad')

# 步骤2：质控
if not ckpt.is_completed('step2_qc'):
    adata_filtered, qc_stats = quality_control(adata)
    adata_filtered.write('./results/02_filtered_data.h5ad')
    ckpt.save_checkpoint('step2_qc',
                        file='./results/02_filtered_data.h5ad')
else:
    print("步骤2已完成，加载数据...")
    adata_filtered = sc.read_h5ad('./results/02_filtered_data.h5ad')

# 获取进度百分比
progress = ckpt.get_progress_percentage()
print(f"完成进度: {progress}%")
```

## 内存监控

断点续传功能与内存监控完全集成：

```bash
# 生成内存使用报告
$ sc-pipeline \
    --sample_info samples.csv \
    --output_dir ./results \
    --resume
```

完成后会生成：
- `memory_usage.csv` - 内存使用详细报告
- `memory_usage.png` - 内存使用趋势图

## 注意事项

### 1. 文件完整性
- 确保中间结果文件没有被删除
- 如果文件损坏，使用 `--reset` 重新运行

### 2. 参数一致性
- 继续运行时应使用相同的参数
- 如果需要修改参数，建议使用 `--reset`

### 3. 断点文件管理
- 断点文件保存在输出目录：`{output_dir}/.checkpoint.json`
- 删除断点文件等同于 `--reset`
- 不同的输出目录有独立的断点

### 4. 磁盘空间
- 中间结果文件会占用磁盘空间
- 完成后可以删除不需要的中间文件
- 建议保留最终的整合数据

## 常见问题

### Q: 断点文件在哪里？
**A:** 在输出目录下的 `.checkpoint.json` 文件

### Q: 如何完全重新运行？
**A:** 使用 `--reset` 参数，或删除 `.checkpoint.json` 文件

### Q: 修改了参数后能继续运行吗？
**A:** 可以，但建议使用 `--reset` 重新运行以保证一致性

### Q: 断点续传会跳过哪些步骤？
**A:** 只跳过已完成的完整步骤，不会跳过部分完成的步骤

### Q: 如果中间文件被删除了怎么办？
**A:** 使用 `--reset` 重新运行，或手动编辑 `.checkpoint.json` 删除对应的步骤记录

### Q: 可以在不同机器上继续运行吗？
**A:** 需要同时复制输出目录和中间结果文件

## 与旧版本的兼容性

如果您使用的是 v1.0 版本，升级到 v2.0 后：
1. 旧的命令行参数完全兼容
2. 新增的断点功能是可选的（不使用 `--resume` 时行为与旧版本相同）
3. 可以安全升级，不影响现有工作流

## 技术实现

断点续传基于以下设计：
1. **CheckpointManager 类**：管理断点状态
2. **JSON 持久化**：断点信息保存为 JSON 文件
3. **步骤标识**：每个主要步骤有唯一标识符
4. **文件验证**：检查中间文件是否存在
5. **异常处理**：捕获错误并保持断点完整性

## 更多信息

- 完整文档: [docs/](docs/)
- 使用示例: [USAGE_EXAMPLES.md](USAGE_EXAMPLES.md)
- 模块化说明: [MODULARIZATION.md](MODULARIZATION.md)
