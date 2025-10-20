#!/bin/bash

# ============================================
# 配置部分
# ============================================
PYTHON_BIN="/File/PATH/TO/YOUR/python"
PY_SCRIPT="/File/PATH/TO/YOUR/sc_pipeline.py"

# 输入输出
SAMPLE_INFO="samples.csv"
OUTPUT_DIR="./results"

# 质控参数
MIN_GENES=200
MAX_GENES=6000
MAX_PCT_MITO=20
DOUBLET_METHOD="scrublet"
DOUBLET_THRESHOLD=0.25

# 整合参数
INTEGRATION_METHOD="harmony"
N_PCS=30
N_NEIGHBORS=10
RESOLUTION=1.1
GENE_NUM=2000
UMAP_MIN_DIST=0.5

# 队列配置
QUEUE="c01"
NUM_CORES=88

# ============================================
# 构建命令
# ============================================
cmd="csub -q ${QUEUE} -R \"span[hosts=1]\" -n ${NUM_CORES} \
    -e sc_pipeline.err -o sc_pipeline.out \
    ${PYTHON_BIN} ${PY_SCRIPT} \
    --sample_info ${SAMPLE_INFO} \
    --output_dir ${OUTPUT_DIR} \
    --min_genes ${MIN_GENES} \
    --max_genes ${MAX_GENES} \
    --max_pct_mito ${MAX_PCT_MITO} \
    --doublet_method ${DOUBLET_METHOD} \
    --doublet_threshold ${DOUBLET_THRESHOLD} \
    --integration_method ${INTEGRATION_METHOD} \
    --n_pcs ${N_PCS} \
    --n_neighbors ${N_NEIGHBORS} \
    --resolution ${RESOLUTION} \
    --gene_num ${GENE_NUM} \
    --umap_min_dist ${UMAP_MIN_DIST}"

# 如果需要细胞类型注释，取消下面的注释
# cmd="${cmd} --run_annotation --celltypist_model /path/to/model.pkl"

# ============================================
# 提交任务
# ============================================
echo "=========================================="
echo "Single Cell RNA-seq Pipeline Job Submission"
echo "=========================================="
echo "Queue: ${QUEUE}"
echo "Cores: ${NUM_CORES}"
echo "Sample info: ${SAMPLE_INFO}"
echo "Output dir: ${OUTPUT_DIR}"
echo "=========================================="
echo ""
echo "Submitting job..."

eval ${cmd}

if [ $? -eq 0 ]; then
    echo ""
    echo "✓ Job submitted successfully!"
    echo ""
    echo "Monitor your job:"
    echo "  - Check queue: bjobs"
    echo "  - View output: tail -f sc_pipeline.out"
    echo "  - View errors: tail -f sc_pipeline.err"
    echo ""
else
    echo ""
    echo "✗ Job submission failed!"
    echo ""
    exit 1
fi
