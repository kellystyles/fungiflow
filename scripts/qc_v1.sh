#!/bin/bash
#SBATCH --array=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --partition=parallel
#SBATCH --time=1-00:00:00
#SBATCH --job-name=qc_v1
#SBATCH -o /nfs/scratch/styleske/fungiflow/scripts/qc_v1_%A_%a.out
#SBATCH -e /nfs/scratch/styleske/fungiflow/scripts/qc_v1_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kelly.styles@vuw.ac.nz

TWD=$(pwd)
qc=/nfs/scratch/styleske/singularity_files/fungiflow/qc_v1.sif

# -------------------------------
# Welcome to the MetaFun workflow
# -------------------------------
#        ~ QC script ~
#
# This script evaluates short (Illumina) reads firstly using FastQC to evaluate reads,
# before parsing the reads for adapter sequences, then trimming the reads using Trimmomatic.
# These trimmed reads are then re-evaluated with FastQC. MultiQC will collate all thee FastQC 
# reports for easy viewing. The FastQC emulator 'Falco' is also available for use in place
# of FastQC. Falco is a bit faster than FastQC so may be the better choice for larger datasets.
# 
# Inputs:
# - raw paired Illumina reads in 'fq.gz' format
# - make sure the variables are set
#
# Please ensure that all raw reads fed into this script are prefixed with a number 
# and underscore like so: "X_" Otherwise, the SLURM task manager won't recognise your files
# To easily do this to a batch of files, use the name_change.sh script in the 'scripts' folder
# 
# Output:
# - adapters can be found in /project_directory/data/adapters/short in 'fasta' format
# - trimmed reads can be found in /project_directory/data/trimmed/short in 'fg.gz' format
# - the FastQC and MultiQC output files can be found in /project_directory/output/(raw_QC and trimmed_QC)/short 

### Variables ###
project="mixed"         # loads project name to target your project directory
eval="fastqc"           # read evaluator (FastQC or Falco)
sif1="qc_v1.sif"        # Singularity image to use for QC
sif2="assembly_v2.sif"  # Singularity image to use for Assembly

###############################
####     Script Begins     ####
###############################

# makes logfile
CD=$(cd .. && pwd)
TWD="$CD"/projects/"$project" # working directory
time=$(date "+%H%M_%d_%m_%Y") # timestamp
log="$time"_"$project"_qc_${SLURM_ARRAY_TASK_ID}.txt
touch "$TWD"/"$log"

# Initialize Singularity module on Raapoi
module load singularity/3.6.1
qc="$CD"/images/"$sif1"
assembly="$CD"/images/"$sif2"
echo "~~~ Quality Control of $project Illumina reads ~~~

The working directory is $TWD
Images:
 - QC       = $qc
 - Assembly = $assembly 
Date: $(date "+%A %d %b %Y")
Time: $(date "+%T")
Software versions:
 - FastQC v0.11.9
 - Falco v0.2.1
 - TrimmomaticPE
 - MultiQC v1.9

Files to be processed:" >> "$TWD"/"$log"
ls -v "$TWD"/data/raw/short/*.fq* > files.txt
while IFS= read -r line; do 
    echo " - ${line##*/}" >> "$TWD"/"$log"
done < files.txt

echo "
     Running QC container
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" >> "$TWD"/"$log"

# ```
# QC
# ```

# fastqc on raw reads
mkdir -p "$TWD"/output/raw/short
(date && echo "| Performing fastqc on raw paired reads") >> "$TWD"/"$log"
singularity exec "$qc" \
"$eval" \
-o "$TWD"/output/raw/short \
"$TWD"/data/raw/short/${SLURM_ARRAY_TASK_ID}_*.fq*

# MultiQC to collect fastqc reports
(date && echo "| Collecting fastqc results for raw paired reads") >> "$TWD"/"$log"
singularity exec "$qc" \
multiqc \
"$TWD"/output/raw/short \
-o "$TWD"/output/raw/short/multiqc

# uses adap_ID.sh to find all adapters in each file and writes them to a file
(date && echo "| Identifying adapter sequences in ${SLURM_ARRAY_TASK_ID}") >> "$TWD"/"$log"
mkdir -p "$TWD"/data/adapters/short
singularity exec "$qc" \
bash /adap_ID.sh \
"$TWD"/data/raw/short/${SLURM_ARRAY_TASK_ID}_R1.fq.gz | grep "^>" -A 1 --no-group-separator > \
"$TWD"/data/adapters/short/${SLURM_ARRAY_TASK_ID}_adapter_file.fasta
singularity exec "$qc" \
bash /adap_ID.sh \
"$TWD"/data/raw/short/${SLURM_ARRAY_TASK_ID}_R2.fq.gz | grep "^>" -A 1 --no-group-separator >> \
"$TWD"/data/adapters/short/${SLURM_ARRAY_TASK_ID}_adapter_file.fasta
    
# appends TruSeq2-PE adapters onto end of created adapter files
singularity exec "$qc" \
cat /TruSeq3-PE-2.fa >> "$TWD"/data/adapters/short/${SLURM_ARRAY_TASK_ID}_adapter_file.fasta

# trims the reads with Trimmomatic
(date && echo "| Trimming Illumina reads no. ${SLURM_ARRAY_TASK_ID} with Trimmomatic") >> "$TWD"/"$log"
mkdir -p "$TWD"/data/trimmed/short
singularity exec "$qc" \
TrimmomaticPE -threads 24 -trimlog "$TWD"/data/raw/short/${SLURM_ARRAY_TASK_ID}_trimmomatic.log \
"$TWD"/data/raw/short/${SLURM_ARRAY_TASK_ID}_R1.fq* "$TWD"/data/raw/short/${SLURM_ARRAY_TASK_ID}_R2.fq* \
"$TWD"/data/trimmed/short/${SLURM_ARRAY_TASK_ID}_forward_trimmed_1P.fq.gz \
"$TWD"/data/trimmed/short/${SLURM_ARRAY_TASK_ID}_forward_trimmed_1U.fq.gz \
"$TWD"/data/trimmed/short/${SLURM_ARRAY_TASK_ID}_reverse_trimmed_2P.fq.gz \
"$TWD"/data/trimmed/short/${SLURM_ARRAY_TASK_ID}_reverse_trimmed_2U.fq.gz \
ILLUMINACLIP:"$TWD"/data/adapters/short/${SLURM_ARRAY_TASK_ID}_adapter_file.fasta:2:30:10:4:4:/true \
TRAILING:9 SLIDINGWINDOW:4:15 MINLEN:36

# fastqc on trimmed reads
mkdir -p "$TWD"/output/trimmed/short
(date && echo "| Performing fastqc on trimmed paired reads") >> "$TWD"/"$log"
singularity exec "$qc" \
"$eval" \
-o "$TWD"/output/trimmed/short \
-a "$TWD"/data/adapters/short/${SLURM_ARRAY_TASK_ID}_adapter_file.fasta \
"$TWD"/data/trimmed/short/${SLURM_ARRAY_TASK_ID}_*P.fq*

# MultiQC to collect fastqc reports
(date && echo "| Collecting fastqc results for trimmed paired reads") >> "$TWD"/"$log"
singularity exec "$qc" \
multiqc \
"$TWD"/output/trimmed/short \
-o "$TWD"/output/trimmed/short/multiqc

(date && echo " ~~ QC Complete ~~") >> "$TWD"/"$log"

echo "
   Running Assembly container
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" >> "$TWD"/"$log"

# ```
# Assembly
# ```

# Assembly with MetaSPADES
(date && echo "| Assembling trimmed reads ${SLURM_ARRAY_TASK_ID} with MetaSPADES ") >> "$TWD"/"$log"
mkdir -p "$TWD"/data/assembly/short/${SLURM_ARRAY_TASK_ID}
singularity exec "$assembly" \
metaspades.py --threads 4 --memory 12 -k 21,33,55,77,99,127 \
--pe1-1 "$TWD"/data/trimmed/short/${SLURM_ARRAY_TASK_ID}_forward_trimmed_1P.fq.gz \
--pe1-2 "$TWD"/data/trimmed/short/${SLURM_ARRAY_TASK_ID}_reverse_trimmed_2P.fq.gz \
-o "$TWD"/data/assembly/short/${SLURM_ARRAY_TASK_ID}
(date && echo "| Finished assembling reads ${SLURM_ARRAY_TASK_ID} with MetaSPADES ") >> "$TWD"/"$log"

echo "
If you found this script useful, or are having issues, send me an email at kelly.styles@vuw.ac.nz
Have a great day :D" >> "$TWD"/"$log"
