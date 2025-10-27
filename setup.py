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
    description="A modular single-cell RNA-seq analysis pipeline with transfer and monitoring capabilities",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/thesecondfox/sc-rna-pipeline",
    project_urls={
        "Bug Tracker": "https://github.com/thesecondfox/sc-rna-pipeline/issues",
        "Documentation": "https://github.com/thesecondfox/sc-rna-pipeline/tree/main/docs",
        "Source Code": "https://github.com/thesecondfox/sc-rna-pipeline",
    },
    packages=find_packages(),
    package_data={
        'sc_pipeline': ['config/*.yaml'],
    },
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
    extras_require={
        'transfer': ['paramiko>=2.7.0'],
        'dev': ['pytest>=6.0', 'black>=21.0', 'flake8>=3.9'],
    },
    entry_points={
        "console_scripts": [
            "sc-pipeline=sc_pipeline.main:main",
        ],
    },
    keywords="single-cell RNA-seq bioinformatics scanpy genomics file-transfer memory-monitoring",
)
