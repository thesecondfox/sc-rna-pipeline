# 安装指南

## 系统要求

### 硬件要求

- **CPU**: 多核处理器（推荐 8 核或以上）
- **内存**: 最少 16GB RAM（推荐 32GB 或更多）
- **存储**: 足够的磁盘空间（取决于数据量，建议 100GB+）

### 软件要求

- **操作系统**: Linux, macOS, 或 Windows (WSL)
- **Python**: 3.8 或更高版本
- **Git**: 用于克隆仓库

## 安装方法

### 方法 1: 使用 Conda（强烈推荐）

Conda 可以更好地管理依赖关系，特别是对于复杂的科学计算包。

#### 步骤 1: 安装 Conda

如果尚未安装 Conda，请选择以下之一：

**Miniconda（推荐，轻量级）**:
```bash
# Linux
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# macOS
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
bash Miniconda3-latest-MacOSX-x86_64.sh
```

#### 步骤 2: 克隆仓库
```bash
git clone https://github.com/thesecondfox/sc-rna-pipeline.git
cd sc-rna-pipeline
```

#### 步骤 3: 创建 Conda 环境
```bash
# 创建新环境
conda create -n sc-pipeline python=3.8

# 激活环境
conda activate sc-pipeline
```

#### 步骤 4: 安装依赖
```bash
# 安装核心依赖
pip install -r requirements.txt

# 或者使用 conda 安装（某些包可能更稳定）
conda install -c conda-forge scanpy python-igraph leidenalg
pip install scrublet harmonypy celltypist
```

### 方法 2: 使用 virtualenv

如果您不使用 Conda，可以使用 Python 的 virtualenv。

#### 步骤 1: 克隆仓库
```bash
git clone https://github.com/thesecondfox/sc-rna-pipeline.git
cd sc-rna-pipeline
```

#### 步骤 2: 创建虚拟环境
```bash
# Linux/macOS
python3 -m venv sc-pipeline-env
source sc-pipeline-env/bin/activate

# Windows
python -m venv sc-pipeline-env
sc-pipeline-env\Scripts\activate
```

#### 步骤 3: 安装依赖
```bash
pip install --upgrade pip
pip install -r requirements.txt
```

## 验证安装

运行以下命令验证安装是否成功：
```bash
python sc_pipeline.py --help
```

如果看到帮助信息，说明安装成功！

您还可以检查关键包的版本：
```bash
python -c "import scanpy; print('scanpy:', scanpy.__version__)"
python -c "import anndata; print('anndata:', anndata.__version__)"
python -c "import scrublet; print('scrublet:', scrublet.__version__)"
```

## 可选组件

### Harmony（批次校正）
```bash
pip install harmonypy
```

### CellTypist（细胞类型注释）
```bash
pip install celltypist

# 下载模型（可选）
celltypist --update-models
```

## 常见问题

### 问题 1: 编译错误

**症状**: 安装时出现 `error: command 'gcc' failed` 或类似错误

**解决方法**:

**Ubuntu/Debian**:
```bash
sudo apt-get update
sudo apt-get install build-essential python3-dev
```

**CentOS/RHEL**:
```bash
sudo yum groupinstall "Development Tools"
sudo yum install python3-devel
```

**macOS**:
```bash
xcode-select --install
```

### 问题 2: igraph 安装失败

**症状**: `ERROR: Failed building wheel for python-igraph`

**解决方法**:

**Ubuntu/Debian**:
```bash
sudo apt-get install libigraph0-dev
pip install python-igraph
```

**macOS (使用 Homebrew)**:
```bash
brew install igraph
pip install python-igraph
```

**使用 Conda**:
```bash
conda install -c conda-forge python-igraph
```

### 问题 3: leidenalg 安装失败

**解决方法**:
```bash
conda install -c conda-forge leidenalg
```

### 问题 4: 依赖冲突

**症状**: 版本冲突错误

**解决方法**:

1. 使用 Conda 而不是 pip
2. 创建新的干净环境
3. 按顺序安装包：
```bash
# 创建新环境
conda create -n sc-pipeline-clean python=3.8
conda activate sc-pipeline-clean

# 按顺序安装
conda install -c conda-forge scanpy python-igraph leidenalg
pip install scrublet harmonypy celltypist
```

### 问题 5: HDF5 错误

**症状**: `HDF5 library version mismatched error`

**解决方法**:
```bash
pip uninstall h5py
pip install --no-cache-dir h5py
```

### 问题 6: macOS Apple Silicon (M1/M2) 问题

**解决方法**:
```bash
# 创建 ARM64 环境
CONDA_SUBDIR=osx-arm64 conda create -n sc-pipeline python=3.8
conda activate sc-pipeline
conda config --env --set subdir osx-arm64
pip install -r requirements.txt
```

## 在 HPC 集群上安装

### 使用模块系统
```bash
# 加载必需的模块
module load python/3.8
module load gcc/9.3.0

# 创建环境
python -m venv sc-pipeline-env
source sc-pipeline-env/bin/activate

# 安装
pip install -r requirements.txt
```

### 使用 Conda 在 HPC
```bash
# 在用户目录安装 Miniconda
cd ~
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3

# 初始化
~/miniconda3/bin/conda init bash
source ~/.bashrc

# 创建环境
conda create -n sc-pipeline python=3.8
conda activate sc-pipeline
pip install -r requirements.txt
```

## 更新

### 更新代码
```bash
cd sc-rna-pipeline
git pull origin main
```

### 更新依赖
```bash
pip install --upgrade -r requirements.txt
```

### 更新 CellTypist 模型
```bash
celltypist --update-models
```

## 卸载

### 删除 Conda 环境
```bash
conda deactivate
conda remove -n sc-pipeline --all
```

### 删除 virtualenv
```bash
deactivate
rm -rf sc-pipeline-env
```

### 删除代码
```bash
rm -rf sc-rna-pipeline
```

## 下一步

安装完成后，请查看：
- [使用说明](usage.md)
- [完整教程](tutorial.md)
- [常见问题](faq.md)

## 获取帮助

如果安装过程中遇到问题：
1. 查看本文档的常见问题部分
2. 搜索 [GitHub Issues](https://github.com/thesecondfox/sc-rna-pipeline/issues)
3. 提交新的 Issue 并附上详细的错误信息
