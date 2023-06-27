#!/bin/bash

# make conda environment
module load Miniconda3/4.9.2
conda create -y -n fungiflow python=3
source activate fungiflow
conda install -y -c conda-forge biopython=1.81 seaborn

# install EggNOG database for Funannotate
mkdir -p /nfs/scratch/campbela/databases/eggnog
export EGGNOG_DATA_DIR=/nfs/scratch/campbel2/databases/eggnog
module load singularity/3.7.3
singularity exec funannotate_v2.sif \
download_eggnog_data.py \
singularity exec funannotate_v2.sif \
create_dbs.py -m diamond --dbname fungi --taxa Fungi