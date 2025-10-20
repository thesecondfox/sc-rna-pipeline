from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="sc-rna-pipeline",
    version="1.0.0",
    author="thesecondfox",
    author_email="thesecondfox@users.noreply.github.com",
    description="A complete single-cell RNA-seq analysis pipeline",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/thesecondfox/sc-rna-pipeline",
    project_urls={
        "Bug Tracker": "https://github.com/thesecondfox/sc-rna-pipeline/issues",
        "Documentation": "https://github.com/thesecondfox/sc-rna-pipeline/tree/main/docs",
        "Source Code": "https://github.com/thesecondfox/sc-rna-pipeline",
    },
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "sc-pipeline=sc_pipeline:main",
        ],
    },
    keywords="single-cell RNA-seq bioinformatics scanpy genomics",
)
```

## 5. .gitignore
```
# Python
__pycache__/
*.py[cod]
*$py.class
*.so
.Python
build/
develop-eggs/
dist/
downloads/
eggs/
.eggs/
lib/
lib64/
parts/
sdist/
var/
wheels/
pip-wheel-metadata/
share/python-wheels/
*.egg-info/
.installed.cfg
*.egg
MANIFEST

# Virtual environments
venv/
ENV/
env/
.venv

# IDE
.vscode/
.idea/
*.swp
*.swo
*~
.DS_Store

# Jupyter Notebook
.ipynb_checkpoints

# Data files (don't commit large data)
*.h5ad
*.csv
*.tsv
*.mtx
*.gz
*.h5
*.loom

# Results
results/
output/
logs/
*.out
*.err

# Temporary files
tmp/
temp/
.cache/

# OS
.DS_Store
.DS_Store?
._*
.Spotlight-V100
.Trashes
ehthumbs.db
Thumbs.db

# Keep example files
!example/*.csv
!example/*.md

# Keep documentation
!docs/**/*

# Keep tests
!tests/**/*
