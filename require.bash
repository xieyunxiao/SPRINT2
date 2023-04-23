#bash

conda config --add channels bioconda
conda config --add channels conda-forge

conda install -c bioconda htslib=1.2.1
conda install -c conda-forge bedtools=2.30.0
conda install -c conda-forge samtools=1.2
conda install -c conda-forge bwa=0.7.17
conda install -c bioconda blat=36x7
