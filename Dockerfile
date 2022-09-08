FROM rocker/r-ver:4.1.2

# Install base utilities
RUN apt-get update && \
    apt-get install -y build-essential && \
    apt-get install -y wget && \
    apt-get install -y libz-dev && \
    apt-get install -y libbz2-dev && \
    apt-get install -y liblzma-dev && \
    apt-get install -y libcurl4-openssl-dev && \
    apt-get install -y libxml2-dev && \
    apt-get install -y libglpk-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# install conda
# ENV CONDA_DIR /opt/conda

# RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
#     /bin/bash ~/miniconda.sh -b -p /opt/conda

# ENV PATH=$CONDA_DIR/bin:$PATH

# RUN conda config --add channels conda-forge
# RUN conda config --add channels bioconda


# Install necessary R packages
RUN Rscript -e 'install.packages("devtools", quiet = TRUE)'
RUN Rscript -e 'install.packages("BiocManager", quiet = TRUE)'
RUN Rscript -e 'BiocManager::install()'
RUN Rscript -e 'install.packages("ArchR")'
RUN Rscript -e 'install.packages("argparse")'


# install scatools
RUN Rscript -e 'devtools::install_github("https://github.com/mjz1/scatools", quiet = TRUE)'