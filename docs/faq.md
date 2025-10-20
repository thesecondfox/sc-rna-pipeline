# 常见问题 (FAQ)

## 安装问题

### Q: 安装时出现编译错误怎么办？

**A**: 确保安装了必要的编译工具：
```bash
# Ubuntu/Debian
sudo apt-get install build-essential python3-dev

# CentOS/RHEL
sudo yum groupinstall "Development Tools"

# macOS
xcode-select --install
```

### Q: python-igraph 安装失败？

**A**: 先安装 igraph C 库：
```bash
# Ubuntu/Debian
sudo apt-get install libigraph0-dev

# macOS
brew install igraph

# 或使用 conda
conda install -c conda-forge python-igraph
```

### Q: 在 HPC 集群上安装遇到权限问题？

**A**: 在用户目录创建虚拟环境：
```bash
cd ~
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
~/miniconda3/bin/conda init bash
source ~/.bashrc

conda create -n sc-pipeline python=3.8
conda activate sc-pipeline
pip install -r requirements.txt
```

## 数据准备问题

### Q: 如何处理 10X Genomics CellRanger 输出？

**A**: CellRanger 输出的 `filtered_feature_bc_matrix` 目录可以直接使用：
```csv
Path,SampleName
/path/to/sample1/outs/filtered_feature_bc_matrix,Sample1
```

### Q: 我的数据是 H5 格式，如何处理？

**A**: 先转换为 MTX 格式：
```python
import scanpy as sc
import scipy.io as io
import pandas as pd

# 读取 H5 文件
adata = sc.read_10x_h5("filtered_feature_bc_matrix.h5")

# 导出为 MTX 格式
io.mmwrite("matrix.mtx", adata.X.T)
pd.DataFrame(adata.obs_names).to_csv("barcodes.tsv", sep="\t", header=False, index=False)
pd.DataFrame(adata.var_names).to_csv("features.tsv", sep="\t", header=False, index=False)
```

### Q: 如何处理多个批次的数据？

**A**: 在样品信息文件中添加 `Batch` 列：
```csv
Path,SampleName,Batch
/data/batch1/sample1,S1,Batch1
/data/batch1/sample2,S2,Batch1
/data/batch2/sample3,S3,Batch2
```

## 运行问题

### Q: 运行时内存不足怎么办？

**A**: 几个解决方案：

1. **减少高变基因数量**:
```bash
--gene_num 1000
```

2. **更严格的质控**:
```bash
--min_genes 500 --max_genes 5000
```

3. **分批处理样品**

4. **使用更大内存的节点**

### Q: 程序运行很慢？

**A**: 

1. **使用更多 CPU 核心**:
```bash
-n 88  # HPC 上增加核心数
```

2. **减少计算量**:
```bash
--gene_num 1000
--n_pcs 20
--doublet_method none
```

### Q: Harmony 整合失败？

**A**: 

1. 检查是否安装了 harmonypy:
```bash
pip install harmonypy
```

2. 尝试使用 Combat:
```bash
--integration_method combat
```

### Q: 出现语法错误？

**A**: 检查 Python 文件末尾，确保没有未闭合的三引号。文件应该以：
```python
if __name__ == "__main__":
    main()
```

## 参数选择问题

### Q: 如何选择合适的质控参数？

**A**: 

1. 先用默认参数运行一次
2. 查看 `qc_statistics.csv`
3. 根据质控统计调整参数

查看分布的示例代码：
```python
import scanpy as sc
import matplotlib.pyplot as plt

adata = sc.read_h5ad("results/raw_data.h5ad")
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

fig, axes = plt.subplots(1, 3, figsize=(15, 4))
axes[0].hist(adata.obs['n_genes_by_counts'], bins=100)
axes[0].set_xlabel('Genes per cell')
axes[1].hist(adata.obs['total_counts'], bins=100)
axes[1].set_xlabel('Counts per cell')
axes[2].hist(adata.obs['pct_counts_mt'], bins=100)
axes[2].set_xlabel('Mitochondrial %')
plt.show()
```

### Q: 如何选择聚类
