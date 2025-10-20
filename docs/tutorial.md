# 完整教程

本教程将通过一个完整的流程，展示如何使用本流程分析单细胞 RNA-seq 数据。

## 数据准备

### 1. 创建项目目录
```bash
mkdir my_project
cd my_project
```

### 2. 准备样品信息文件

创建 `samples.csv`:
```csv
Path,SampleName,Group,Condition
/path/to/sample1,Sample1,Control,Untreated
/path/to/sample2,Sample2,Treatment,Drug_A
/path/to/sample3,Sample3,Control,Untreated
/path/to/sample4,Sample4,Treatment,Drug_A
```

### 3. 验证数据完整性
```bash
# 检查每个样品目录
for path in /path/to/sample*; do
    echo "Checking $path"
    ls -lh $path
done
```

## 运行分析

### 基础分析
```bash
python /path/to/sc_pipeline.py \
    --sample_info samples.csv \
    --output_dir ./results \
    --min_genes 200 \
    --max_genes 6000 \
    --max_pct_mito 20 \
    --integration_method harmony
```

### HPC 任务提交

创建并编辑 `submit.sh`:
```bash
#!/bin/bash

PYTHON_BIN="python"
PY_SCRIPT="/path/to/sc_pipeline.py"

cmd="csub -q c01 -n 88 \
    -e analysis.err -o analysis.out \
    ${PYTHON_BIN} ${PY_SCRIPT} \
    --sample_info samples.csv \
    --output_dir ./results \
    --integration_method harmony"

eval ${cmd}
```

提交任务：
```bash
chmod +x submit.sh
./submit.sh
```

## 结果查看

### 1. 查看输出文件
```bash
ls -lh results/
```

### 2. 查看质控统计
```bash
cat results/qc_statistics.csv
```

### 3. 查看 UMAP 图

如果在服务器上，需要下载图片到本地查看：
```bash
scp user@server:~/my_project/results/*.png ./local_dir/
```

## 下游分析

### 1. 加载数据
```python
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

# 读取整合后的数据
adata = sc.read_h5ad("results/integrated_data.h5ad")

print(f"细胞数: {adata.n_obs:,}")
print(f"基因数 (.X): {adata.n_vars:,}")
print(f"基因数 (.raw): {adata.raw.n_vars:,}")
```

### 2. 基本可视化
```python
# UMAP 可视化
sc.pl.umap(adata, color='leiden')
sc.pl.umap(adata, color=['Group', 'Condition'])

# 查看 QC 指标
sc.pl.violin(adata, ['n_genes', 'n_counts', 'pct_counts_mt'], 
             groupby='leiden', rotation=90)
```

### 3. 寻找 Marker 基因
```python
# 计算每个 cluster 的 marker 基因
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

# 可视化
sc.pl.rank_genes_groups(adata, n_genes=20)

# 保存结果
markers = sc.get.rank_genes_groups_df(adata, group=None)
markers.to_csv("results/marker_genes.csv", index=False)
```

### 4. 可视化特定基因
```python
# 定义感兴趣的基因
genes_of_interest = ['CD3D', 'CD8A', 'CD4', 'CD19', 'MS4A1']

# 在 UMAP 上显示
sc.pl.umap(adata, color=genes_of_interest, ncols=3)

# 小提琴图
sc.pl.violin(adata, genes_of_interest, groupby='leiden')

# 点图
sc.pl.dotplot(adata, genes_of_interest, groupby='leiden')
```

### 5. 细胞类型注释
```python
# 基于 marker 基因手动注释
cluster_annotation = {
    '0': 'T cells',
    '1': 'B cells',
    '2': 'Monocytes',
    '3': 'NK cells',
    # 根据实际情况添加...
}

adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_annotation)
adata.obs['cell_type'] = adata.obs['cell_type'].fillna('Unknown')

# 可视化注释结果
sc.pl.umap(adata, color='cell_type', legend_loc='on data')

# 保存注释后的数据
adata.write("results/annotated_data.h5ad")
```

### 6. 差异表达分析
```python
# 比较两个组之间的差异
sc.tl.rank_genes_groups(adata, 'Group', method='wilcoxon',
                        groups=['Treatment'], reference='Control')

# 获取差异基因
deg = sc.get.rank_genes_groups_df(adata, group='Treatment')
deg_sig = deg[deg['pvals_adj'] < 0.05]

print(f"发现 {len(deg_sig)} 个显著差异基因")

# 保存结果
deg_sig.to_csv("results/DEGs_Treatment_vs_Control.csv", index=False)

# 可视化
sc.pl.rank_genes_groups_violin(adata, groups='Treatment', n_genes=10)
```

### 7. 组成分析
```python
# 统计每个组的细胞类型组成
comp = pd.crosstab(adata.obs['cell_type'], adata.obs['Group'], 
                   normalize='columns') * 100

# 可视化
comp.plot(kind='bar', stacked=True, figsize=(10, 6))
plt.ylabel('Percentage (%)')
plt.xlabel('Cell Type')
plt.title('Cell Type Composition')
plt.legend(title='Group')
plt.tight_layout()
plt.savefig('results/composition_analysis.png', dpi=300)
```

## 高级分析

### 1. 轨迹分析
```python
# 使用 PAGA
sc.tl.paga(adata, groups='cell_type')
sc.pl.paga(adata, save='_paga.png')

# 计算伪时间
sc.tl.diffmap(adata)
sc.tl.dpt(adata)
sc.pl.umap(adata, color='dpt_pseudotime')
```

### 2. 子集分析
```python
# 提取特定细胞类型进行深入分析
t_cells = adata[adata.obs['cell_type'] == 'T cells'].copy()

# 重新聚类
sc.tl.leiden(t_cells, resolution=0.5)
sc.pl.umap(t_cells, color='leiden')

# 寻找亚群 marker
sc.tl.rank_genes_groups(t_cells, 'leiden')
```

## 生成报告
```python
# 创建分析摘要
summary = {
    'Total cells': adata.n_obs,
    'Total genes': adata.raw.n_vars,
    'Number of clusters': adata.obs['leiden'].nunique(),
    'Samples': adata.obs['SampleName'].unique().tolist(),
}

with open('results/analysis_summary.txt', 'w') as f:
    for key, value in summary.items():
        f.write(f"{key}: {value}\n")
```

## 下一步

1. **功能富集分析**: 使用 GO/KEGG
2. **细胞-细胞通讯**: 使用 CellPhoneDB
3. **整合多组学**: 结合 ATAC-seq
4. **发表准备**: 优化图表质量

## 参考资料

- [Scanpy Documentation](https://scanpy.readthedocs.io/)
- [Single-cell Best Practices](https://www.sc-best-practices.org/)
- [Seurat Tutorials](https://satijalab.org/seurat/)
