# 使用示例

本文档展示如何使用模块化后的单细胞RNA-seq分析流程。

## 1. 基本命令行使用

### 1.1 基础分析流程

```bash
sc-pipeline \
    --sample_info samples.csv \
    --output_dir ./results \
    --min_genes 200 \
    --max_genes 6000 \
    --max_pct_mito 20 \
    --integration_method harmony
```

### 1.2 启用内存监控

```bash
sc-pipeline \
    --sample_info samples.csv \
    --output_dir ./results \
    --enable_memory_monitor \
    --memory_threshold 80.0
```

### 1.3 包含细胞类型注释

```bash
sc-pipeline \
    --sample_info samples.csv \
    --output_dir ./results \
    --run_annotation \
    --celltypist_model /path/to/model.pkl
```

## 2. 模块化使用（Python脚本）

### 2.1 基本流程

```python
from sc_pipeline import (
    read_sample,
    quality_control,
    data_integration,
    save_h5ad
)

# 读取数据
adata = read_sample('samples.csv')

# 质控
adata_filtered, qc_stats = quality_control(
    adata,
    min_genes=200,
    max_genes=6000,
    max_pct_mito=20,
    return_stats=True
)

# 数据整合
adata_integrated = data_integration(
    adata_filtered,
    method='harmony',
    n_pcs=30,
    resolution=1.1,
    output_dir='./results'
)

# 保存结果
save_h5ad(adata_integrated, './results/integrated_data.h5ad')
```

### 2.2 使用内存监控

```python
from sc_pipeline import MemoryMonitor, memory_profiler, log_memory

# 方法1: 使用装饰器
@memory_profiler
def analyze_data():
    adata = read_sample('samples.csv')
    adata_filtered, _ = quality_control(adata)
    return adata_filtered

# 方法2: 手动监控
monitor = MemoryMonitor(
    threshold_percent=80.0,
    check_interval=2.0
)

# 启动后台监控
monitor.start_monitoring()

# 执行分析
adata = read_sample('samples.csv')

# 记录内存使用
log_memory("读取数据后: ")

# 质控
adata_filtered, _ = quality_control(adata)

log_memory("质控后: ")

# 停止监控并打印摘要
monitor.stop_monitoring()
monitor.print_summary()

# 方法3: 设置阈值回调
def memory_warning(mem_info):
    print(f"⚠️ 内存使用率过高: {mem_info['system_percent']:.1f}%")
    print(f"可用内存: {mem_info['system_available_mb']:.1f} MB")

monitor = MemoryMonitor(threshold_percent=75.0)
monitor.set_threshold_callback(memory_warning)
monitor.start_monitoring()

# ... 执行分析 ...

monitor.stop_monitoring()
```

### 2.3 文件传输功能

#### 2.3.1 本地文件操作

```python
from sc_pipeline import copy_file, calculate_checksum, get_file_info

# 拷贝文件（带进度）
def progress_callback(transferred, total):
    percent = transferred / total * 100
    print(f"\r进度: {percent:.1f}%", end='')

copy_file(
    src='/path/to/large_file.h5ad',
    dst='/backup/large_file.h5ad',
    progress_callback=progress_callback,
    verify=True
)

# 计算文件校验和
checksum = calculate_checksum('/path/to/file.h5ad')
print(f"MD5: {checksum}")

# 获取文件信息
info = get_file_info('/path/to/file.h5ad')
print(f"文件大小: {info['size_mb']:.2f} MB")
```

#### 2.3.2 远程文件传输（SFTP）

```python
from sc_pipeline import SFTPTransfer

# 创建SFTP连接
sftp = SFTPTransfer(
    host='remote.server.com',
    port=22,
    username='your_username',
    key_file='/path/to/private_key'  # 或使用 password='your_password'
)

# 连接服务器
sftp.connect()

# 上传文件
sftp.upload(
    local_path='./results/integrated_data.h5ad',
    remote_path='/remote/path/integrated_data.h5ad'
)

# 下载文件
sftp.download(
    remote_path='/remote/path/data.h5ad',
    local_path='./data/data.h5ad'
)

# 列出远程目录
files = sftp.list_dir('/remote/path')
print(f"远程文件: {files}")

# 断开连接
sftp.disconnect()
```

#### 2.3.3 断点续传

```python
from sc_pipeline import TransferManager, SFTPTransfer

# 创建传输管理器
manager = TransferManager(
    chunk_size=1024 * 1024 * 10,  # 10MB分块
    state_dir='./.transfer_state',
    max_retries=3,
    retry_delay=5
)

# 创建SFTP传输对象
sftp = SFTPTransfer(
    host='remote.server.com',
    port=22,
    username='username',
    password='password'
)
sftp.connect()

# 定义传输函数
def upload_chunk(src, dst, offset, length):
    """上传文件块"""
    sftp.upload(src, f"{dst}.part{offset}")

# 执行断点续传上传
success = manager.upload_file(
    local_path='./large_file.h5ad',
    remote_path='/remote/large_file.h5ad',
    transfer_func=upload_chunk,
    resume=True,  # 启用断点续传
    verify=True
)

if success:
    print("文件上传成功！")
else:
    print("文件上传失败，可以稍后重试")

# 查看传输任务状态
tasks = manager.list_tasks()
for task in tasks:
    print(f"任务: {task['src']} -> {task['dst']}")
    print(f"状态: {task['status']}")
    print(f"进度: {task['transferred_bytes']} / {task['total_bytes']}")

# 清理已完成的任务
manager.clean_completed_tasks()

sftp.disconnect()
```

### 2.4 仅使用特定模块

```python
# 仅使用IO模块
from sc_pipeline.io import readbgi, save_h5ad

adata = readbgi('/path/to/data', 'sample1')
save_h5ad(adata, 'sample1.h5ad')

# 仅使用预处理模块
from sc_pipeline.preprocessing import quality_control

adata_filtered, stats = quality_control(adata, return_stats=True)

# 仅使用整合模块
from sc_pipeline.integration import data_integration

adata_integrated = data_integration(
    adata_filtered,
    method='harmony',
    output_dir='./results'
)

# 仅使用注释模块
from sc_pipeline.annotation import celltypist_annotation

adata_annotated = celltypist_annotation(
    adata_integrated,
    model_path='/path/to/model.pkl'
)
```

## 3. 高级用法

### 3.1 自定义分析流程

```python
from sc_pipeline import (
    read_sample,
    quality_control,
    data_integration,
    MemoryMonitor,
    save_h5ad,
    save_csv
)

# 启动内存监控
monitor = MemoryMonitor(threshold_percent=80.0)
monitor.start_monitoring()

# 读取数据
print("步骤1: 读取数据...")
adata = read_sample('samples.csv')
print(f"读取 {adata.n_obs} 个细胞, {adata.n_vars} 个基因")

# 质控
print("\n步骤2: 质量控制...")
adata_filtered, qc_stats = quality_control(
    adata,
    min_genes=200,
    max_genes=6000,
    max_pct_mito=20,
    doublet_method='scrublet',
    return_stats=True
)

# 保存质控统计
save_csv(qc_stats, './results/qc_stats.csv', index=False)

# 数据整合
print("\n步骤3: 数据整合...")
adata_integrated = data_integration(
    adata_filtered,
    method='harmony',
    n_pcs=30,
    resolution=1.1,
    output_dir='./results'
)

# 保存结果
print("\n保存结果...")
save_h5ad(adata_integrated, './results/final_data.h5ad')

# 停止监控并显示摘要
monitor.stop_monitoring()
print("\n内存使用摘要:")
monitor.print_summary()
```

### 3.2 批量处理多个数据集

```python
import os
from sc_pipeline import (
    read_sample,
    quality_control,
    data_integration,
    save_h5ad,
    MemoryMonitor
)

# 数据集列表
datasets = [
    'dataset1/samples.csv',
    'dataset2/samples.csv',
    'dataset3/samples.csv'
]

monitor = MemoryMonitor()
monitor.start_monitoring()

for i, dataset in enumerate(datasets):
    print(f"\n处理数据集 {i+1}/{len(datasets)}: {dataset}")

    # 创建输出目录
    output_dir = f'./results/dataset_{i+1}'
    os.makedirs(output_dir, exist_ok=True)

    # 读取和处理
    adata = read_sample(dataset)
    adata_filtered, _ = quality_control(adata)
    adata_integrated = data_integration(
        adata_filtered,
        output_dir=output_dir
    )

    # 保存
    save_h5ad(adata_integrated, f'{output_dir}/integrated.h5ad')

    # 清理内存
    del adata, adata_filtered, adata_integrated
    import gc
    gc.collect()

monitor.stop_monitoring()
monitor.print_summary()
```

## 4. 配置文件使用

可以创建配置文件 `config.yaml` 来管理参数：

```yaml
# config.yaml
input:
  sample_info: samples.csv

output:
  directory: ./results

qc:
  min_genes: 200
  max_genes: 6000
  max_pct_mito: 20.0
  doublet_method: scrublet
  doublet_threshold: 0.25

integration:
  method: harmony
  n_pcs: 30
  n_neighbors: 10
  resolution: 1.1
  gene_num: 2000
  umap_min_dist: 0.5

monitoring:
  enable: true
  threshold: 80.0
```

然后在Python中加载：

```python
import yaml
from sc_pipeline import *

# 加载配置
with open('config.yaml', 'r') as f:
    config = yaml.safe_load(f)

# 使用配置
adata = read_sample(config['input']['sample_info'])

adata_filtered, _ = quality_control(
    adata,
    min_genes=config['qc']['min_genes'],
    max_genes=config['qc']['max_genes'],
    max_pct_mito=config['qc']['max_pct_mito']
)

adata_integrated = data_integration(
    adata_filtered,
    method=config['integration']['method'],
    n_pcs=config['integration']['n_pcs'],
    output_dir=config['output']['directory']
)
```

## 5. 故障排除

### 5.1 内存不足

```python
# 使用内存监控识别问题
from sc_pipeline import MemoryMonitor

monitor = MemoryMonitor(threshold_percent=70.0)

def high_memory_warning(mem_info):
    print(f"⚠️ 内存使用率: {mem_info['system_percent']:.1f}%")
    if mem_info['system_percent'] > 80:
        print("建议: 减少处理的样品数量或使用更大的内存")

monitor.set_threshold_callback(high_memory_warning)
monitor.start_monitoring()

# 执行分析...
```

### 5.2 传输中断

```python
# 断点续传会自动处理中断
# 只需重新运行相同的命令即可继续

from sc_pipeline import TransferManager

manager = TransferManager()

# 第一次尝试（可能中断）
try:
    manager.upload_file(local_path='large.h5ad',
                       remote_path='/remote/large.h5ad',
                       transfer_func=upload_func,
                       resume=True)
except:
    print("传输中断")

# 第二次尝试（从断点继续）
manager.upload_file(local_path='large.h5ad',
                   remote_path='/remote/large.h5ad',
                   transfer_func=upload_func,
                   resume=True)  # 会自动从上次断点继续
```

## 6. 更多信息

- 完整文档: [docs/](docs/)
- 安装指南: [docs/installation.md](docs/installation.md)
- 教程: [docs/tutorial.md](docs/tutorial.md)
- FAQ: [docs/faq.md](docs/faq.md)
