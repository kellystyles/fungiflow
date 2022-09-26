#!/bin/bash
#SBATCH --array=403,331,213,216,212,217,402,330,187,182,338,172,177,173,186,183,316,313,159,234,231,235,230,317,158,312,238,155,154,151,239,258,259,254,373,376,377,250,255,358,359,273,351,286,355,350,277,272,241,244,363,366,362,367,240,245,248,266,296,341,292,297,345,348,349,197,192,162,328,167,60,329,163,166,196,193,321,324,203,202,207,320,325,228,43,229,306,149,303,224,221,225,220,307,302,299,268,298,260,265,342,347,295,343,346,294,261,264,369,247,242,365,360,364,361,243,300,305,222,227,301,304,309,308,327,322,168,205,200,198,204,201,199,326,169,191,194,209,164,161,165,63,190,195,208,153,319,156,318,152,157,310,232,237,233,311,314,58,219,181,184,72,174,171,73,175,218,185,337,178,188,215,210,189,214,211,336,333,179,401,270,275,280,352,357,353,356,271,274,278,288,252,380,375,370,374,371,256,253,379,378%6
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --partition=parallel
#SBATCH --time=8:00:00
#SBATCH --job-name=fungiflow
#SBATCH -o /nfs/scratch/styleske/fungiflow/fungiflow_%A_%a.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kelly.styles@vuw.ac.nz

module load old-mod-system/Miniconda3/4.9.2
module load singularity/3.7.3

source activate fungiflow

python3 /nfs/scratch/styleske/fungiflow/fungiflow.py -d "/nfs/scratch/styleske/fungiflow/projects/microcera_test" -a ${SLURM_ARRAY_TASK_ID} \
-if ${SLURM_ARRAY_TASK_ID}_1.fq.gz \
-ir ${SLURM_ARRAY_TASK_ID}_2.fq.gz \
-c ${SLURM_CPUS_PER_TASK} -m ${SLURM_MEM_PER_NODE} -ant -t "isolate" -b \
-s /nfs/scratch/styleske/singularity_files/fungiflow_v0-4_latest.sif \
-sf /nfs/scratch/styleske/singularity_files/funannotate_gm.sif \
-sa /nfs/scratch/styleske/singularity_files/antismash_6.0.1--pyhdfd78af_0.sif \
-idb /nfs/scratch/styleske/databases/fungi_ITS/ITS_RefSeq_Fungi \
-kdb /nfs/scratch/styleske/databases/kraken2_std/ \
-bdb /nfs/scratch/styleske/databases/NCBI_nt/nt
