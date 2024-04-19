FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="f07de3c28593a9890bfed8e9d79692f967c8ab07be59bb34323c5de216119b7f"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: env/bash.yml
#   prefix: /conda-envs/6d084b213928bc5754913e32cd1935bc
#   channels:
#     - conda-forge
#   dependencies:
#     - bash=5.1.16
#     - pigz=2.6
#     - bc=1.07.1
RUN mkdir -p /conda-envs/6d084b213928bc5754913e32cd1935bc
COPY env/bash.yml /conda-envs/6d084b213928bc5754913e32cd1935bc/environment.yaml

# Conda environment:
#   source: env/bedtools.yml
#   prefix: /conda-envs/4451a440ce51beb4e0d4fe61d6f741cb
#   channels:
#     - bioconda
#   dependencies:
#     - bedtools=2.30.0
RUN mkdir -p /conda-envs/4451a440ce51beb4e0d4fe61d6f741cb
COPY env/bedtools.yml /conda-envs/4451a440ce51beb4e0d4fe61d6f741cb/environment.yaml

# Conda environment:
#   source: env/bowtie2-samtools.yml
#   prefix: /conda-envs/d6da1c6a6f0dcffab87c2f07c238cd21
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - bowtie2=2.5.0
#     - samtools=1.16.1
#     - pigz=2.6
RUN mkdir -p /conda-envs/d6da1c6a6f0dcffab87c2f07c238cd21
COPY env/bowtie2-samtools.yml /conda-envs/d6da1c6a6f0dcffab87c2f07c238cd21/environment.yaml

# Conda environment:
#   source: env/csaw.yml
#   prefix: /conda-envs/1421d9bbfeb09cc410b4b1f430b7df16
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - bioconductor-csaw=1.32.0
#     - bioconductor-summarizedexperiment=1.28.0
RUN mkdir -p /conda-envs/1421d9bbfeb09cc410b4b1f430b7df16
COPY env/csaw.yml /conda-envs/1421d9bbfeb09cc410b4b1f430b7df16/environment.yaml

# Conda environment:
#   source: env/deeptools.yml
#   prefix: /conda-envs/d657bc7dbf87a49dc8eb27c608db8a30
#   channels:
#     - bioconda
#   dependencies:
#     - deeptools=3.5.1
RUN mkdir -p /conda-envs/d657bc7dbf87a49dc8eb27c608db8a30
COPY env/deeptools.yml /conda-envs/d657bc7dbf87a49dc8eb27c608db8a30/environment.yaml

# Conda environment:
#   source: env/macs2.yml
#   prefix: /conda-envs/ca1d466f3efa2d6f714b374ba60a8bae
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - macs2=2.2.7.1
RUN mkdir -p /conda-envs/ca1d466f3efa2d6f714b374ba60a8bae
COPY env/macs2.yml /conda-envs/ca1d466f3efa2d6f714b374ba60a8bae/environment.yaml

# Conda environment:
#   source: env/tidyverse.yml
#   prefix: /conda-envs/74938091e56ba9db67bd24adfec8cd1b
#   channels:
#     - r
#     - conda-forge
#   dependencies:
#     - r-tidyverse
RUN mkdir -p /conda-envs/74938091e56ba9db67bd24adfec8cd1b
COPY env/tidyverse.yml /conda-envs/74938091e56ba9db67bd24adfec8cd1b/environment.yaml

# Conda environment:
#   source: env/wget.yml
#   prefix: /conda-envs/5283609becff4ef3b65849f54fe8f2c2
#   channels:
#     - anaconda
#   dependencies:
#     - wget=1.21.3
RUN mkdir -p /conda-envs/5283609becff4ef3b65849f54fe8f2c2
COPY env/wget.yml /conda-envs/5283609becff4ef3b65849f54fe8f2c2/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v1.20.0/bio/bowtie2/align/environment.yaml
#   prefix: /conda-envs/a0d814b6b7bab88fa10bae3c06e51492
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - bowtie2 =2.5.0
#     - samtools =1.16.1
#     - snakemake-wrapper-utils =0.5.0
RUN mkdir -p /conda-envs/a0d814b6b7bab88fa10bae3c06e51492
ADD https://github.com/snakemake/snakemake-wrappers/raw/v1.20.0/bio/bowtie2/align/environment.yaml /conda-envs/a0d814b6b7bab88fa10bae3c06e51492/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v1.20.0/bio/bowtie2/build/environment.yaml
#   prefix: /conda-envs/be2fcc5a63c2af8ef0faff1a6ad17cbb
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - bowtie2 =2.5.0
RUN mkdir -p /conda-envs/be2fcc5a63c2af8ef0faff1a6ad17cbb
ADD https://github.com/snakemake/snakemake-wrappers/raw/v1.20.0/bio/bowtie2/build/environment.yaml /conda-envs/be2fcc5a63c2af8ef0faff1a6ad17cbb/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v1.20.0/bio/deeptools/computematrix/environment.yaml
#   prefix: /conda-envs/03fe262236eea767b423afbe31e6a282
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - deeptools =3.5.1
RUN mkdir -p /conda-envs/03fe262236eea767b423afbe31e6a282
ADD https://github.com/snakemake/snakemake-wrappers/raw/v1.20.0/bio/deeptools/computematrix/environment.yaml /conda-envs/03fe262236eea767b423afbe31e6a282/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v1.20.0/bio/fastqc/environment.yaml
#   prefix: /conda-envs/8678da093d0504bc9ac36dfd48d4c62e
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - fastqc =0.11.9
RUN mkdir -p /conda-envs/8678da093d0504bc9ac36dfd48d4c62e
ADD https://github.com/snakemake/snakemake-wrappers/raw/v1.20.0/bio/fastqc/environment.yaml /conda-envs/8678da093d0504bc9ac36dfd48d4c62e/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v1.20.0/bio/multiqc/environment.yaml
#   prefix: /conda-envs/eebc5adeb5e51053b4346bc959a40041
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - multiqc =1.13
RUN mkdir -p /conda-envs/eebc5adeb5e51053b4346bc959a40041
ADD https://github.com/snakemake/snakemake-wrappers/raw/v1.20.0/bio/multiqc/environment.yaml /conda-envs/eebc5adeb5e51053b4346bc959a40041/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v1.20.0/bio/picard/markduplicates/environment.yaml
#   prefix: /conda-envs/e96c80318da18d792a199b6b90d2db3a
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - picard =2.27.4
#     - samtools =1.16.1
#     - snakemake-wrapper-utils =0.5.0
RUN mkdir -p /conda-envs/e96c80318da18d792a199b6b90d2db3a
ADD https://github.com/snakemake/snakemake-wrappers/raw/v1.20.0/bio/picard/markduplicates/environment.yaml /conda-envs/e96c80318da18d792a199b6b90d2db3a/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v1.20.0/bio/samtools/idxstats/environment.yaml
#   prefix: /conda-envs/bb3906be70e50c58a2f9c73826e6c973
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - samtools =1.16.1
#     - snakemake-wrapper-utils =0.5.0
RUN mkdir -p /conda-envs/bb3906be70e50c58a2f9c73826e6c973
ADD https://github.com/snakemake/snakemake-wrappers/raw/v1.20.0/bio/samtools/idxstats/environment.yaml /conda-envs/bb3906be70e50c58a2f9c73826e6c973/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v1.20.0/bio/samtools/index/environment.yaml
#   prefix: /conda-envs/851a5df29edaf302e173deb7719ea8ba
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - samtools =1.16.1
RUN mkdir -p /conda-envs/851a5df29edaf302e173deb7719ea8ba
ADD https://github.com/snakemake/snakemake-wrappers/raw/v1.20.0/bio/samtools/index/environment.yaml /conda-envs/851a5df29edaf302e173deb7719ea8ba/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v1.20.0/bio/trimmomatic/se/environment.yaml
#   prefix: /conda-envs/1ca4ed0b0eff4b7ba35d2c223e15abfe
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - trimmomatic ==0.36
#     - pigz ==2.3.4
#     - snakemake-wrapper-utils ==0.1.3
RUN mkdir -p /conda-envs/1ca4ed0b0eff4b7ba35d2c223e15abfe
ADD https://github.com/snakemake/snakemake-wrappers/raw/v1.20.0/bio/trimmomatic/se/environment.yaml /conda-envs/1ca4ed0b0eff4b7ba35d2c223e15abfe/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v1.21.1/bio/homer/annotatePeaks/environment.yaml
#   prefix: /conda-envs/e0275b642c24c0c7be8fd493ac926809
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - homer =4.11
RUN mkdir -p /conda-envs/e0275b642c24c0c7be8fd493ac926809
ADD https://github.com/snakemake/snakemake-wrappers/raw/v1.21.1/bio/homer/annotatePeaks/environment.yaml /conda-envs/e0275b642c24c0c7be8fd493ac926809/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v1.21.2/bio/bamtools/filter_json/environment.yaml
#   prefix: /conda-envs/fade5d588c8d41f538d33a69e79736fb
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - bamtools =2.5.2
RUN mkdir -p /conda-envs/fade5d588c8d41f538d33a69e79736fb
ADD https://github.com/snakemake/snakemake-wrappers/raw/v1.21.2/bio/bamtools/filter_json/environment.yaml /conda-envs/fade5d588c8d41f538d33a69e79736fb/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/6d084b213928bc5754913e32cd1935bc --file /conda-envs/6d084b213928bc5754913e32cd1935bc/environment.yaml && \
    mamba env create --prefix /conda-envs/4451a440ce51beb4e0d4fe61d6f741cb --file /conda-envs/4451a440ce51beb4e0d4fe61d6f741cb/environment.yaml && \
    mamba env create --prefix /conda-envs/d6da1c6a6f0dcffab87c2f07c238cd21 --file /conda-envs/d6da1c6a6f0dcffab87c2f07c238cd21/environment.yaml && \
    mamba env create --prefix /conda-envs/1421d9bbfeb09cc410b4b1f430b7df16 --file /conda-envs/1421d9bbfeb09cc410b4b1f430b7df16/environment.yaml && \
    mamba env create --prefix /conda-envs/d657bc7dbf87a49dc8eb27c608db8a30 --file /conda-envs/d657bc7dbf87a49dc8eb27c608db8a30/environment.yaml && \
    mamba env create --prefix /conda-envs/ca1d466f3efa2d6f714b374ba60a8bae --file /conda-envs/ca1d466f3efa2d6f714b374ba60a8bae/environment.yaml && \
    mamba env create --prefix /conda-envs/74938091e56ba9db67bd24adfec8cd1b --file /conda-envs/74938091e56ba9db67bd24adfec8cd1b/environment.yaml && \
    mamba env create --prefix /conda-envs/5283609becff4ef3b65849f54fe8f2c2 --file /conda-envs/5283609becff4ef3b65849f54fe8f2c2/environment.yaml && \
    mamba env create --prefix /conda-envs/a0d814b6b7bab88fa10bae3c06e51492 --file /conda-envs/a0d814b6b7bab88fa10bae3c06e51492/environment.yaml && \
    mamba env create --prefix /conda-envs/be2fcc5a63c2af8ef0faff1a6ad17cbb --file /conda-envs/be2fcc5a63c2af8ef0faff1a6ad17cbb/environment.yaml && \
    mamba env create --prefix /conda-envs/03fe262236eea767b423afbe31e6a282 --file /conda-envs/03fe262236eea767b423afbe31e6a282/environment.yaml && \
    mamba env create --prefix /conda-envs/8678da093d0504bc9ac36dfd48d4c62e --file /conda-envs/8678da093d0504bc9ac36dfd48d4c62e/environment.yaml && \
    mamba env create --prefix /conda-envs/eebc5adeb5e51053b4346bc959a40041 --file /conda-envs/eebc5adeb5e51053b4346bc959a40041/environment.yaml && \
    mamba env create --prefix /conda-envs/e96c80318da18d792a199b6b90d2db3a --file /conda-envs/e96c80318da18d792a199b6b90d2db3a/environment.yaml && \
    mamba env create --prefix /conda-envs/bb3906be70e50c58a2f9c73826e6c973 --file /conda-envs/bb3906be70e50c58a2f9c73826e6c973/environment.yaml && \
    mamba env create --prefix /conda-envs/851a5df29edaf302e173deb7719ea8ba --file /conda-envs/851a5df29edaf302e173deb7719ea8ba/environment.yaml && \
    mamba env create --prefix /conda-envs/1ca4ed0b0eff4b7ba35d2c223e15abfe --file /conda-envs/1ca4ed0b0eff4b7ba35d2c223e15abfe/environment.yaml && \
    mamba env create --prefix /conda-envs/e0275b642c24c0c7be8fd493ac926809 --file /conda-envs/e0275b642c24c0c7be8fd493ac926809/environment.yaml && \
    mamba env create --prefix /conda-envs/fade5d588c8d41f538d33a69e79736fb --file /conda-envs/fade5d588c8d41f538d33a69e79736fb/environment.yaml && \
    mamba clean --all -y
