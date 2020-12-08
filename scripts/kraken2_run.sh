#!/bin/bash
#SBATCH --array=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --partition=parallel
#SBATCH --time=1-00:00:00
#SBATCH --job-name=kraken2_taxonomy
#SBATCH -o /nfs/scratch/styleske/fungiflow/scripts/kraken2_%A_%a.out
#SBATCH -e /nfs/scratch/styleske/fungiflow/scripts/kraken2_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kelly.styles@vuw.ac.nz

# This script evaluates files for OTUs using Kraken2

# adds directory info for files to be evaluated
project="mixed"
# adds directory info for type of files to be evaluated
type="raw" # "raw", "trimmed", "assembled"
# adds length info for reads
length="short"
# are the reads paired?
paired="yes"
# adds kraken2 database to use
db="kraken2_standard"
# Singularity image to use
sif="kraken2_v1.sif"

###############################
####     Script Begins     ####
###############################

# creates log file
CD=$(cd .. && pwd)
TWD="$CD"/projects/"$project" # working directory
time=$(date "+%H%M_%d_%m_%Y") # timestamp
log="$time"_kraken2_"$project".txt
touch "$TWD"/"$log"

init=$(date +%s)

module load singularity/3.6.1

# begin Kraken2 analyses

ver="v2.1.0"

echo "~~~ Kraken2 Analyses of $project $type data ~~~

The working directory is $TWD
Date: $(date "+%A %d %b %Y")
Time: $(date "+%T")
Software versions:nano krake    
 - Kraken2 $ver

Files to be processed:" >>"$TWD"/"$log"
for f in "$TWD"/data/"$type"/"$length"/*; do
    echo " - ${f##*/}" >>"$TWD"/"$log"
done
echo "" >>"$TWD"/"$log"

   
# Kraken 2 analysis
(date "+%T| Beggining Kraken2 analysis of file for array ${SLURM_ARRAY_TASK_ID}") >> "$TWD"/"$log"
mkdir -p "$TWD"/data/classified/
cat <<eof > $TWD/data/classified/execute_${SLURM_ARRAY_TASK_ID}.sh
#!/bin/bash
singularity exec $CD/images/$sif \\
kraken2 \\
--db $CD/databases/$db \\
--threads 4 \\
--use-names \\
eof

case $paired in 
    "yes") echo "--paired \\
    --classified-out "$TWD"/data/classified/${SLURM_ARRAY_TASK_ID}_cseqs#.fq \\
    --unclassified-out "$TWD"/data/classified/${SLURM_ARRAY_TASK_ID}_useqs#.fq \\
    "$TWD"/data/"$type"/short/${SLURM_ARRAY_TASK_ID}_R1*.fq* \\
    "$TWD"/data/"$type"/short/${SLURM_ARRAY_TASK_ID}_R2*.fq*" ;;    
    *) echo "--classified-out "$TWD"/data/classified/${SLURM_ARRAY_TASK_ID}_cseqs.fq \\
    --unclassified-out "$TWD"/data/classified/${SLURM_ARRAY_TASK_ID}_useqs.fq \\
    "$TWD"/data/raw/"$length"/${SLURM_ARRAY_TASK_ID}_*.f*" ;;
esac >> "$TWD"/data/classified/execute_${SLURM_ARRAY_TASK_ID}.sh

source "$TWD"/data/classified/execute_${SLURM_ARRAY_TASK_ID}.sh

comp=$(date +%s)
time=$((init-comp))
hours=$((time / 3600))
minutes=$(( (time % 3600) / 60 ))
seconds=$(( (time % 3600) % 60 ))

# compiling kraken2 outputs into one file

(date "+%T| Compiling Kraken2 output into a single file") >> "$TWD"/"$log"
touch "$TWD"/classified/"$project"_kraken2_output_data.txt
for file in "$TWD"/classified/*.fq; do
    echo " | Procesing output of $file" >> "$TWD"/"$log"
    name=${file##*/}
    awk '/kraken/ {print $0, $name}' >>  "$TWD"/classified/"$project"_kraken2_output_data.txt
done
(date "+%T| Finished processing output data") >> "$TWD"/"$log"

echo "
Total runtime: $hours:$minutes:$seconds (hh:mm:ss)

If you found this script useful, or are having issues, send me an email at kelly.styles@vuw.ac.nz
Have a great day :D" >> "$TWD"/"$log"
