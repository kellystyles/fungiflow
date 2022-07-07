#!/bin/bash

# Installing shell packages
apt-get update
apt-get -y install \
build-essential manpages-dev apt-utils locales wget default-jre liblzma-dev \
libboost-all-dev perl wget build-essential python3 python3-pip \
zlib1g-dev pkg-config libfreetype6-dev libpng-dev python-matplotlib \
openjdk-11-jdk openjdk-11-jre trimmomatic
apt-get clean

# Installing Miniconda
if [ ! -d /usr/local/anaconda ]; then
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    -O ~/conda.sh && bash ~/conda.sh -bfp /usr/local
    rm ~/conda.sh
fi
conda update conda -y -q
. /usr/local/etc/profile.d/conda.sh
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels defaults

# Install micromaba
conda create -n environ mamba -c conda-forge -c defaults
conda activate environ
mamba install -y -c conda-forge -c bioconda -y kraken2 python=3.8 pilon spades flye masurca fmlrc \
racon medaka krona itsx fastqc nanocomp minimap2 seaborn docopt tqdm wget pyyaml git pysam blobtools
    
mamba update -y kraken2 pilon spades flye masurca fmlrc \
racon medaka krona itsx fastqc nanocomp minimap2 seaborn docopt tqdm wget pyyaml git blobtools
mamba install -y -c bioconda pysam --update-deps
mamba clean -a

# Installing specific software
# - BBMap suite
wget https://versaweb.dl.sourceforge.net/project/bbmap/BBMap_38.96.tar.gz 
tar xvzf BBMap_38.96.tar.gz && rm BBMap_38.96.tar.gz
mv bbmap/* /usr/bin/

# - NCBI+ suite
VERSION=2.11.0
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/$VERSION/ncbi-blast-$VERSION+-x64-linux.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/$VERSION/ncbi-blast-$VERSION+-x64-linux.tar.gz.md5
md5sum -c ncbi-blast-2.11.0+-x64-linux.tar.gz.md5
tar zxvpf ncbi-blast-$VERSION+-x64-linux.tar.gz && rm ncbi-blast-$VERSION+-x64-linux.tar.gz
mv ncbi-blast-2.11.0+/bin/* /usr/bin/

# - Quast
wget https://downloads.sourceforge.net/project/quast/quast-5.0.2.tar.gz
tar -xzf quast-5.0.2.tar.gz
cd quast-5.0.2
./setup.py install
QUAST_JSON="/usr/local/envs/environ/lib/python3.8/site-packages/quast-5.0.2-py3.8.egg/quast_libs/site_packages/jsontemplate/jsontemplate.py"
sed -i 's/import cgi/import html/g' "$QUAST_JSON"
sed -i 's/cgi.escape/html.escape/g' "$QUAST_JSON"

# - BlobToolKit
dir=$()
cd /usr/local/envs/environ/lib/python3.*/site-packages/blobtools/
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz* -P data/
tar zxf data/taxdump.tar.gz -C data/ nodes.dmp names.dmp
blobtools nodesdb --nodes data/nodes.dmp --names data/names.dmp

# Running test scriptlet to confirm installations
python3 env_test_scriptlet.py