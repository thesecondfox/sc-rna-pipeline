# 模块化重构说明

## 概述

本项目已完成模块化重构，将原来的单体文件 `sc_pipeline.py` 拆分为多个功能模块，并新增了文件传输和内存监控功能。

## 新的项目结构

```
sc-rna-pipeline/
├── sc_pipeline/                    # 主包目录
│   ├── __init__.py                # 包初始化文件
│   ├── main.py                    # 命令行入口
│   │
│   ├── io/                        # IO模块
│   │   ├── __init__.py
│   │   ├── reader.py             # 数据读取 (readbgi, read_sample)
│   │   └── writer.py             # 数据写入 (save_h5ad, save_csv)
│   │
│   ├── preprocessing/             # 预处理模块
│   │   ├── __init__.py
│   │   └── qc.py                 # 质量控制 (quality_control)
│   │
│   ├── integration/               # 整合模块
│   │   ├── __init__.py
│   │   └── integration.py        # 数据整合 (data_integration)
│   │
│   ├── annotation/                # 注释模块
│   │   ├── __init__.py
│   │   └── celltype.py           # 细胞类型注释 (celltypist_annotation)
│   │
│   ├── transfer/                  # 文件传输模块 (新增)
│   │   ├── __init__.py
│   │   ├── local.py              # 本地文件操作
│   │   ├── remote.py             # 远程文件传输 (SFTP/FTP)
│   │   └── manager.py            # 传输管理器 (支持断点续传)
│   │
│   └── utils/                     # 工具模块 (新增)
│       ├── __init__.py
│       └── memory_monitor.py     # 内存监控
│
├── setup.py                       # 安装配置（已更新）
├── requirements.txt               # 依赖列表（已更新）
├── USAGE_EXAMPLES.md             # 使用示例（新增）
├── MODULARIZATION.md             # 本文档
└── docs/                         # 文档目录
```

## 新增功能

### 1. 内存监控

提供实时内存监控功能，帮助用户跟踪和管理内存使用。

**主要特性：**
- 实时监控进程和系统内存使用
- 设置内存阈值告警
- 记录内存使用历史
- 提供装饰器和手动监控两种方式

**使用方法：**

```python
from sc_pipeline import MemoryMonitor, memory_profiler

# 方法1: 装饰器
@memory_profiler
def my_analysis():
    # 分析代码
    pass

# 方法2: 手动监控
monitor = MemoryMonitor(threshold_percent=80.0)
monitor.start_monitoring()
# ... 执行分析 ...
monitor.stop_monitoring()
monitor.print_summary()
```

**命令行使用：**

```bash
sc-pipeline \
    --sample_info samples.csv \
    --enable_memory_monitor \
    --memory_threshold 80.0
```

### 2. 文件传输（支持断点续传）

提供本地和远程文件传输功能，支持断点续传。

**主要特性：**
- 本地文件拷贝（带进度显示）
- 远程文件传输（SFTP/FTP）
- 断点续传（分块传输，支持中断恢复）
- 传输状态持久化
- 自动重试机制
- 文件完整性验证

**使用方法：**

```python
from sc_pipeline import TransferManager, SFTPTransfer

# 创建传输管理器
manager = TransferManager(
    chunk_size=1024*1024*10,  # 10MB分块
    max_retries=3,
    retry_delay=5
)

# 创建SFTP连接
sftp = SFTPTransfer(
    host='remote.server.com',
    username='user',
    key_file='/path/to/key'
)
sftp.connect()

# 上传文件（支持断点续传）
def upload_func(src, dst, offset, length):
    sftp.upload(src, dst)

success = manager.upload_file(
    local_path='large_file.h5ad',
    remote_path='/remote/large_file.h5ad',
    transfer_func=upload_func,
    resume=True  # 启用断点续传
)

# 查看传输状态
tasks = manager.list_tasks()
for task in tasks:
    print(f"{task['src']} -> {task['dst']}: {task['status']}")

sftp.disconnect()
```

## 模块化优势

### 1. 可维护性
- 代码按功能分组，易于理解和维护
- 每个模块职责单一，降低耦合度
- 便于单元测试

### 2. 可扩展性
- 新增功能时只需添加新模块
- 不影响现有模块的稳定性
- 支持插件式开发

### 3. 可复用性
- 各模块可独立使用
- 方便集成到其他项目
- 支持按需导入

### 4. 灵活性
- 用户可以选择性使用模块
- 支持自定义工作流
- 便于与现有代码集成

## 向后兼容性

原有的命令行接口保持不变：

```bash
# 旧版使用方式仍然有效
sc-pipeline --sample_info samples.csv --output_dir ./results
```

原有的 `sc_pipeline.py` 文件保留，可以继续使用。新增功能通过模块化方式提供。

## 迁移指南

### 从单体文件迁移到模块化

**原代码：**
```python
from sc_pipeline import readbgi, quality_control, data_integration
```

**新代码：**
```python
# 方式1: 从主包导入（推荐）
from sc_pipeline import readbgi, quality_control, data_integration

# 方式2: 从子模块导入
from sc_pipeline.io import readbgi
from sc_pipeline.preprocessing import quality_control
from sc_pipeline.integration import data_integration
```

两种方式都兼容，推荐使用方式1（从主包导入）。

## 新增依赖

- `psutil>=5.8.0` - 系统监控（必需）
- `paramiko>=2.7.0` - SFTP传输（可选，安装：`pip install sc-rna-pipeline[transfer]`）

## 性能改进

模块化后的性能特点：

1. **内存管理**：通过内存监控，可以及时发现和处理内存问题
2. **可中断性**：支持断点续传，大文件传输不怕中断
3. **按需加载**：只导入需要的模块，减少启动时间

## 测试

运行测试（如果安装了开发依赖）：

```bash
# 安装开发依赖
pip install -e ".[dev]"

# 运行测试
pytest tests/

# 代码格式检查
flake8 sc_pipeline/

# 代码格式化
black sc_pipeline/
```

## 文档

- [安装指南](docs/installation.md)
- [使用教程](docs/tutorial.md)
- [使用示例](USAGE_EXAMPLES.md)
- [FAQ](docs/faq.md)

## 更新日志

### v1.0.0 - 模块化重构

**新增功能：**
- ✨ 模块化架构
- ✨ 内存监控功能
- ✨ 文件传输模块（支持断点续传）
- ✨ SFTP/FTP远程传输
- ✨ 传输管理器（状态持久化、自动重试）

**改进：**
- 🔧 重构代码结构
- 📝 完善文档和使用示例
- 🧪 添加类型提示
- 🎨 改进代码可读性

**向后兼容：**
- ✅ 保持原有命令行接口
- ✅ 保持原有函数签名
- ✅ 保留原有 sc_pipeline.py

## 贡献

欢迎贡献代码！请遵循以下步骤：

1. Fork 本仓库
2. 创建特性分支 (`git checkout -b feature/AmazingFeature`)
3. 提交更改 (`git commit -m 'Add some AmazingFeature'`)
4. 推送到分支 (`git push origin feature/AmazingFeature`)
5. 创建 Pull Request

## 许可证

MIT License

## 联系方式

如有问题或建议，请提交 Issue。
