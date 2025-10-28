# 断点续传功能 - 快速开始

## 🎯 核心特性

增强版单细胞RNA-seq分析流程提供了强大的**断点续传**功能，让您的分析更可靠、更高效：

✅ **自动保存进度** - 每个步骤完成后自动保存
✅ **智能恢复** - 从中断点无缝继续运行
✅ **数据验证** - 自动检查文件完整性
✅ **错误追踪** - 详细的错误日志和诊断
✅ **交互式操作** - 友好的用户界面和提示

---

## 🚀 快速开始

### 1. 首次运行

```bash
python sc_pipeline_enhanced.py \
    --sample_info samples.csv \
    --output_dir ./results \
    --min_genes 500 \
    --integration_method harmony
```

### 2. 如果中断，从断点恢复

```bash
python sc_pipeline_enhanced.py \
    --sample_info samples.csv \
    --output_dir ./results \
    --resume
```

系统会显示：
```
检测到之前的进度:
================================================================================
✅ 步骤1: 数据读取 (已完成)
✅ 步骤2: 质量控制 (已完成)
❌ 步骤3: 数据整合 (失败)
⏳ 步骤4: 细胞注释 (待执行)
================================================================================

是否从上次断点继续运行？(y/n):
```

### 3. 查看当前进度

```bash
python sc_pipeline_enhanced.py \
    --sample_info samples.csv \
    --output_dir ./results \
    --show-progress
```

---

## 📋 命令行参数

### 断点控制

| 参数 | 说明 | 示例 |
|------|------|------|
| `--resume` | 交互式恢复 | `--resume` |
| `--auto-resume` | 自动恢复（无需确认） | `--auto-resume` |
| `--show-progress` | 只显示进度，不运行 | `--show-progress` |
| `--reset` | 重置所有断点 | `--reset` |

### 分析参数

完整参数列表参见 `--help`：

```bash
python sc_pipeline_enhanced.py --help
```

---

## 📊 工作流程

```
┌─────────────────────────────────────────────────────────────┐
│                     开始分析                                │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ▼
         ┌───────────────────────┐
         │   检查是否有断点？     │
         └─────┬─────────┬───────┘
               │ 有      │ 无
               ▼         ▼
         ┌─────────┐   ┌──────────────┐
         │显示进度  │   │ 从头开始执行  │
         │提示恢复  │   │              │
         └────┬────┘   └──────┬───────┘
              │               │
              └───────┬───────┘
                      ▼
            ┌──────────────────┐
            │  步骤1: 数据读取  │──┐
            └─────────┬────────┘  │
                      ▼            │ 每步完成
            ┌──────────────────┐  │ 自动保存
            │  步骤2: 质量控制  │  │ 断点
            └─────────┬────────┘  │
                      ▼            │
            ┌──────────────────┐  │
            │  步骤3: 数据整合  │◄─┘
            └─────────┬────────┘
                      ▼
            ┌──────────────────┐
            │  步骤4: 细胞注释  │
            └─────────┬────────┘
                      ▼
            ┌──────────────────┐
            │   分析完成！      │
            └──────────────────┘
```

---

## 💡 使用场景

### 场景 1: 内存不足

**问题**: 在数据整合步骤内存溢出

```bash
# 调整参数后继续
python sc_pipeline_enhanced.py \
    --sample_info samples.csv \
    --output_dir ./results \
    --gene_num 2000 \
    --n_pcs 30 \
    --resume
```

### 场景 2: 参数调优

**问题**: 想尝试不同的聚类分辨率

```bash
# 删除步骤3的输出，重新运行
rm ./results/03_integrated_data.h5ad

python sc_pipeline_enhanced.py \
    --sample_info samples.csv \
    --output_dir ./results \
    --resolution 1.2 \
    --resume
```

### 场景 3: 批处理脚本

**问题**: 需要在服务器上无人值守运行

```bash
#!/bin/bash
# auto_analysis.sh

python sc_pipeline_enhanced.py \
    --sample_info samples.csv \
    --output_dir ./results \
    --auto-resume \
    --integration_method harmony

# 如果失败，使用降级参数重试
if [ $? -ne 0 ]; then
    echo "使用降级参数重试..."
    python sc_pipeline_enhanced.py \
        --sample_info samples.csv \
        --output_dir ./results \
        --auto-resume \
        --gene_num 2000
fi
```

---

## 📁 生成的文件

### 断点文件

```
./results/
├── .checkpoint.json       # 断点状态文件
├── .error_log.txt        # 错误日志（如有）
├── 01_raw_data.h5ad      # 步骤1输出
├── 02_filtered_data.h5ad # 步骤2输出
├── 03_integrated_data.h5ad # 步骤3输出
└── 04_annotated_data.h5ad # 步骤4输出（可选）
```

### 断点文件示例

`.checkpoint.json`:
```json
{
  "step1_read": {
    "completed": true,
    "timestamp": "2025-10-28T10:30:15",
    "valid": true,
    "file": "./results/01_raw_data.h5ad",
    "n_cells": 45230
  }
}
```

---

## 🔍 进度查询

### 命令行查询

```bash
python sc_pipeline_enhanced.py \
    --sample_info samples.csv \
    --output_dir ./results \
    --show-progress
```

### Python脚本查询

```python
from sc_pipeline_enhanced import CheckpointManager

# 创建管理器
ckpt = CheckpointManager("./results")

# 显示进度
ckpt.print_status()

# 获取最后完成的步骤
last_step = ckpt.get_last_completed_step()
print(f"最后完成: {last_step}")

# 获取下一步
next_step = ckpt.get_next_step()
print(f"下一步: {next_step}")
```

---

## 🧪 测试

运行功能测试：

```bash
# 交互式测试
python test_checkpoint.py

# 选择测试：
# 1. 基本断点保存和恢复
# 2. 失败后断点恢复
# 3. 断点完整性验证
# 4. 断点重置
```

---

## 📖 详细文档

- **[CHECKPOINT_GUIDE.md](CHECKPOINT_GUIDE.md)** - 完整使用指南
- **[CHANGELOG_CHECKPOINT.md](CHANGELOG_CHECKPOINT.md)** - 详细更新日志
- **[test_checkpoint.py](test_checkpoint.py)** - 测试脚本

---

## 🆚 版本对比

| 功能 | 基础版 | 增强版 |
|------|--------|--------|
| 断点保存 | ✅ | ✅ |
| 文件完整性检查 | ❌ | ✅ |
| 交互式恢复 | ❌ | ✅ |
| 自动恢复 | ❌ | ✅ |
| 错误日志 | ❌ | ✅ |
| 进度查询 | ❌ | ✅ |
| 失败追踪 | ❌ | ✅ |
| 友好提示 | ⚠️ | ✅ |

---

## ⚡ 性能优化

### 1. 内存管理

查看内存使用报告：
```bash
cat ./results/memory_usage.csv
open ./results/memory_usage.png
```

### 2. 参数调优

如遇内存问题，降低这些参数：
- `--gene_num`: 3000 → 2000
- `--n_pcs`: 50 → 30
- `--resolution`: 0.8 → 0.6

### 3. 批处理优化

```bash
# 使用所有CPU核心
python sc_pipeline_enhanced.py --n_jobs -1 ...

# 或指定核心数
python sc_pipeline_enhanced.py --n_jobs 8 ...
```

---

## 🐛 故障排查

### 问题 1: 断点文件损坏

```bash
# 删除断点文件，但保留数据
rm ./results/.checkpoint.json

# 系统会自动重建断点
python sc_pipeline_enhanced.py --resume ...
```

### 问题 2: 数据文件丢失

系统会自动检测并标记为无效：
```
⚠️ 警告: 步骤 step2_qc 的数据文件不存在
```

解决方案：该步骤会自动重新执行

### 问题 3: 查看错误详情

```bash
# 查看错误日志
cat ./results/.error_log.txt

# 查看断点状态
cat ./results/.checkpoint.json | jq
```

---

## 💬 常见问题

**Q: 可以手动修改断点吗？**

A: 可以。编辑 `.checkpoint.json` 文件，修改或删除特定步骤的断点。

**Q: 如何强制重新运行某个步骤？**

A: 删除该步骤的输出文件和断点记录：
```bash
rm ./results/02_filtered_data.h5ad
# 编辑 .checkpoint.json，删除 step2_qc 条目
```

**Q: 断点文件可以跨版本使用吗？**

A: 增强版向后兼容基础版的断点文件。

**Q: 支持并行运行多个分析吗？**

A: 支持。使用不同的 `--output_dir` 即可：
```bash
python sc_pipeline_enhanced.py --output_dir ./results_v1 ...
python sc_pipeline_enhanced.py --output_dir ./results_v2 ...
```

---

## 📞 获取帮助

```bash
# 查看所有参数
python sc_pipeline_enhanced.py --help

# 查看当前版本
python sc_pipeline_enhanced.py --version
```

---

## 📝 示例工作流

### 完整示例

```bash
# 1. 首次运行
python sc_pipeline_enhanced.py \
    --sample_info samples.csv \
    --output_dir ./results \
    --min_genes 500 \
    --max_genes 8000 \
    --integration_method harmony \
    --run_annotation

# 2. 假设在步骤3失败...查看进度
python sc_pipeline_enhanced.py \
    --sample_info samples.csv \
    --output_dir ./results \
    --show-progress

# 3. 查看错误日志
cat ./results/.error_log.txt

# 4. 修复问题后继续
python sc_pipeline_enhanced.py \
    --sample_info samples.csv \
    --output_dir ./results \
    --resume

# 5. 分析完成！查看结果
ls -lh ./results/
```

---

## 🎉 总结

增强版断点续传功能让您：

- 💪 **更可靠** - 不怕中断，随时恢复
- ⚡ **更高效** - 无需重复计算
- 🔍 **更透明** - 详细的进度和错误信息
- 😊 **更友好** - 清晰的提示和交互

**立即开始使用！**

```bash
python sc_pipeline_enhanced.py \
    --sample_info your_samples.csv \
    --output_dir ./your_results
```

---

**版本**: 2.2
**最后更新**: 2025-10-28
**作者**: Claude Code
