# 断点续传功能更新日志

## 版本 2.2 - 增强版断点续传 (2025-10-28)

### 主要改进

#### 1. **增强的CheckpointManager类**

新增功能：
- ✅ **数据完整性验证**: 自动检查保存的数据文件是否存在和完整
- ✅ **文件大小校验**: 记录并验证文件大小，检测损坏的数据
- ✅ **智能恢复判断**: 自动识别有效和无效的断点
- ✅ **详细进度显示**: 展示每个步骤的完成状态、时间和关键统计

#### 2. **完善的错误处理机制**

新增功能：
- ✅ **异常捕获和记录**: 捕获所有错误并保存详细追踪信息
- ✅ **错误日志文件**: 独立的 `.error_log.txt` 记录所有错误
- ✅ **失败步骤标记**: 在断点文件中标记失败的步骤
- ✅ **安全执行包装器**: `safe_step_execution()` 函数统一处理错误

#### 3. **交互式恢复功能**

新增功能：
- ✅ **智能恢复提示**: 检测到断点时自动提示用户
- ✅ **进度预览**: 恢复前显示详细的当前进度
- ✅ **用户确认机制**: 交互式确认是否从断点继续
- ✅ **自动恢复模式**: `--auto-resume` 参数支持无人值守运行

#### 4. **新增命令行参数**

```bash
--resume           # 从断点恢复（交互式）
--auto-resume      # 自动从断点恢复（无需确认）
--reset            # 重置所有断点，从头开始
--show-progress    # 查看当前进度（不运行）
```

### 技术实现

#### CheckpointManager 核心方法

```python
# 数据完整性验证
def _validate_checkpoints(checkpoints):
    """检查文件存在性和大小"""

# 增强的断点保存
def save_checkpoint(step, **kwargs):
    """自动记录文件大小和时间戳"""

# 失败标记
def mark_step_failed(step, error):
    """记录错误信息到断点和日志文件"""

# 智能恢复判断
def should_resume():
    """交互式或自动判断是否恢复"""

# 进度查询
def get_last_completed_step():
    """获取最后完成的步骤"""

def get_next_step():
    """获取下一个待执行步骤"""
```

#### 安全执行包装器

```python
def safe_step_execution(checkpoint_manager, step_name, func, *args, **kwargs):
    """
    统一的步骤执行包装器

    功能：
    - 捕获所有异常（包括KeyboardInterrupt）
    - 自动记录错误到断点管理器
    - 打印详细的错误追踪信息
    - 友好的错误提示
    """
```

### 断点文件格式

#### .checkpoint.json 结构

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
    "error_message": "Unable to allocate array",
    "traceback": "Traceback (most recent call last):\n..."
  }
}
```

#### .error_log.txt 格式

```
================================================================================
步骤: step3_integrate
时间: 2025-10-28T10:50:18.456789
错误类型: MemoryError
错误信息: Unable to allocate array with shape (50000, 3000) and data type float64
详细追踪:
Traceback (most recent call last):
  File "sc_pipeline_enhanced.py", line 1250, in integration_step
    sc.pp.scale(adata, max_value=10)
  ...
MemoryError: Unable to allocate array with shape (50000, 3000) and data type float64
```

### 使用示例

#### 示例 1: 正常运行并中断

```bash
# 首次运行
python sc_pipeline_enhanced.py \
    --sample_info samples.csv \
    --output_dir ./results

# 假设在步骤2失败...

# 修复问题后继续
python sc_pipeline_enhanced.py \
    --sample_info samples.csv \
    --output_dir ./results \
    --resume
```

输出:
```
检测到之前的进度:
================================================================================
当前分析进度
================================================================================
✅ 步骤1: 数据读取
   完成时间: 2025-10-28 10:30:15
   文件: 01_raw_data.h5ad (1250.5 MB)
   细胞数: 45,230

❌ 步骤2: 质量控制 - 失败
   错误: Invalid sample path
   时间: 2025-10-28 10:35:20

⏳ 步骤3: 数据整合 - 待执行
⏳ 步骤4: 细胞注释 - 待执行
================================================================================

是否从上次断点继续运行？(y/n): y

🔄 从断点恢复，跳过已完成的步骤...
```

#### 示例 2: 自动恢复（用于脚本）

```bash
#!/bin/bash
# 自动化分析脚本

python sc_pipeline_enhanced.py \
    --sample_info samples.csv \
    --output_dir ./results \
    --auto-resume \
    --min_genes 500 \
    --integration_method harmony

# 如果失败，调整参数后重试
if [ $? -ne 0 ]; then
    echo "尝试降低参数重新运行..."
    python sc_pipeline_enhanced.py \
        --sample_info samples.csv \
        --output_dir ./results \
        --auto-resume \
        --gene_num 2000 \
        --n_pcs 30
fi
```

#### 示例 3: 查看进度

```bash
python sc_pipeline_enhanced.py \
    --sample_info samples.csv \
    --output_dir ./results \
    --show-progress
```

输出:
```
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

✅ 步骤3: 数据整合
   完成时间: 2025-10-28 11:05:18
   文件: 03_integrated_data.h5ad (850.7 MB)
   细胞数: 38,567
   聚类数: 15

⏳ 步骤4: 细胞注释 - 待执行
================================================================================

📍 下一步: 步骤4: 细胞注释
```

### 优势对比

| 功能 | 原版本 | 增强版本 |
|------|--------|---------|
| 断点保存 | ✅ | ✅ |
| 文件完整性检查 | ❌ | ✅ |
| 错误详细记录 | ❌ | ✅ |
| 交互式恢复 | ❌ | ✅ |
| 自动恢复模式 | ❌ | ✅ |
| 进度查看 | ❌ | ✅ |
| 失败步骤标记 | ❌ | ✅ |
| 错误日志文件 | ❌ | ✅ |
| 用户友好提示 | ⚠️ | ✅ |

### 错误恢复场景

#### 场景 1: 内存不足

**问题**: 步骤3因内存不足失败

**检测**:
```
❌ 步骤3: 数据整合 - 失败
   错误: MemoryError: Unable to allocate array
```

**解决**:
```bash
# 降低参数后继续
python sc_pipeline_enhanced.py \
    --sample_info samples.csv \
    --output_dir ./results \
    --gene_num 2000 \
    --n_pcs 30 \
    --resume
```

#### 场景 2: 数据文件损坏

**问题**: 断点文件存在但数据文件被删除

**检测**:
```
⚠️  警告: 步骤 step2_qc 的数据文件不存在: ./results/02_filtered_data.h5ad
```

**解决**: 系统自动标记为无效，重新执行该步骤

#### 场景 3: 参数调整

**问题**: 需要调整质控参数

**解决**:
```bash
# 删除步骤2的输出
rm ./results/02_filtered_data.h5ad

# 手动编辑断点文件，或直接使用 --reset 重新开始指定步骤
python sc_pipeline_enhanced.py \
    --sample_info samples.csv \
    --output_dir ./results \
    --min_genes 300 \
    --max_pct_mito 20 \
    --resume
```

### 测试

提供了完整的测试脚本 `test_checkpoint.py`:

```bash
# 运行所有测试
python test_checkpoint.py

# 测试内容：
# 1. 基本断点保存和恢复
# 2. 失败后断点恢复
# 3. 断点完整性验证
# 4. 断点重置
```

### 文档

新增文档文件：
- `CHECKPOINT_GUIDE.md`: 详细的使用指南
- `CHANGELOG_CHECKPOINT.md`: 本更新日志
- `test_checkpoint.py`: 功能测试脚本

### 向后兼容性

✅ **完全兼容**: 增强版完全兼容原有的命令行参数和工作流程
- 不使用新参数时，行为与原版本一致
- 原有的断点文件格式仍然支持
- 可以无缝升级

### 未来计划

可能的进一步改进：
- [ ] 支持并行步骤的断点管理
- [ ] 添加断点的版本控制
- [ ] 支持远程断点同步
- [ ] 提供Web界面查看进度
- [ ] 自动性能优化建议

### 反馈和贡献

如有问题或建议，请提交Issue或Pull Request。

---

**版本**: 2.2
**发布日期**: 2025-10-28
**作者**: Claude Code
**许可**: MIT
