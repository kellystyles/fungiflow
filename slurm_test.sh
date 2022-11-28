#!/bin/bash
#SBATCH --array=149,151,152,155,156,157,158,164,166,167,168,173,178,182,186,188,189,190,191,192,195,197,201,202,204,205,210,212,213,215,216,217,219,220,221,228,230,231,233,235,237,238,240,241,248,250,253,256,264,265,266,270,271,273,274,275,280,288,297,298,299,301,303,307,308,312,313,324,325,326,329,330,341,347,348,349,351,352,353,356,357,359,360,361,364,365,367,369,373,374,377,380,60,63,72,73%6
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

python3 /nfs/scratch/styleske/fungiflow/fungiflow.py \
-d "/nfs/scratch/styleske/endolichenic_fungi/all" \
-a ${SLURM_ARRAY_TASK_ID} \
-if ${SLURM_ARRAY_TASK_ID}_1.fq.gz \
-ir ${SLURM_ARRAY_TASK_ID}_2.fq.gz \
-c ${SLURM_CPUS_PER_TASK} -m ${SLURM_MEM_PER_NODE} -ant -t "isolate" \
-s /nfs/scratch/styleske/singularity_files/fungiflow_v0-4_latest.sif \
-sf /nfs/scratch/styleske/singularity_files/funannotate_v0.1.sif \
-sa /nfs/scratch/styleske/singularity_files/antismash_6.0.1--pyhdfd78af_0.sif \
-data /nfs/scratch/styleske/databases/ \
-its -k -b -e