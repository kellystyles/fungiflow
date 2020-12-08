#!/bin/bash -e
#SBATCH --array=16-20%4
#SBATCH --cpus-per-task=24
#SBATCH --mem=320G
#SBATCH --constraint=Intel
#SBATCH --partition=bigmem
#SBATCH --time=1-00:00:00
#SBATCH --job-name=evaluations
#SBATCH -o /nfs/scratch/styleske/fungiflow/scripts/evaluations_%A_%a.out
#SBATCH -e /nfs/scratch/styleske/fungiflow/scripts/evaluations_%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kelly.styles@vuw.ac.nz

### Variables ###
project="lichen_metagenomes"        # loads project name to target your project directory
type="metagenomic"                  # what type of sample data are you assembling ('isolate' or a 'metagenomic')
cpus="24"                           # Number of cpus to use per task
mem="80"                            # Amount of memory per task in GB
eval="fastqc"                       # read evaluator ('fastqc' or 'falco')
paired="yes"                        # are the reads paired?
stddb="kraken2_standard"            # the Kraken2 standard database
fundb="kraken2_fungi"               # the Kraken2 fungi database
kronadb="krona"                     # the Krona taxonomy database name
sif1="qc_v1.sif"                    # Singularity image to use for QC
sif2="kraken2_v1.sif"               # Singularity image to use for Kraken2
sif3="assembly_v2.sif"              # Singularity image to use for Assembly
sif4="antismash_5_1_ver3"           # Singularity image to use of antiSMASH
sif5="itsx_v1.sif"                  # Singularity image to use for ITSx
space_saver="on"                    # removes intermediate files ('on' or 'off')

###############################
####     Script Begins     ####
###############################

start=$(date +%s)

# makes logfile
CD=$(cd .. && pwd)
TWD="$CD"/projects/"$project" # working directory
time=$(date "+%H%M_%d_%m_%Y") # timestamp
log="$time"_"$project"_evaluations_${SLURM_ARRAY_TASK_ID}.txt
touch "$TWD"/"$log"

# Initialize Singularity module on Raapoi
module purge
module load singularity/3.6.1
qc="$CD"/images/"$sif1"
kraken2="$CD"/images/"$sif2"
assembly="$CD"/images/"$sif3"
antismash="$CD"/images/"$sif4"
itsx="$CD"/images/"$sif5"

echo "~~~ Analysis of $project ${SLURM_ARRAY_TASK_ID} Illumina reads ~~~

The working directory is: $TWD

Job name:    $SLURM_JOB_NAME - $project
Job ID:      $SLURM_JOB_ID
Job array:   $SLURM_ARRAY_TASK_ID 

     Assembly Evaluations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" >> "$TWD"/"$log"

final_genome="$TWD"/data/assembly/${SLURM_ARRAY_TASK_ID}/${SLURM_ARRAY_TASK_ID}_scaffold.fasta

# Kraken2 analysis of contigs
(date "+%T| Beggining Kraken2 analysis of contigs for array ${SLURM_ARRAY_TASK_ID}") >> "$TWD"/"$log"
singularity exec "$kraken2" \
kraken2 \
--db "$CD"/databases/"$fundb" \
--threads "$cpus" \
--use-names \
--output "$TWD"/output/classified/${SLURM_ARRAY_TASK_ID}_kraken2_contigs_output.txt \
"$final_genome"

# Generating Kronogram of Kraken2 output
(date "+%T| Generating Kronogram of Kraken2 contig analysis for array ${SLURM_ARRAY_TASK_ID}") >> "$TWD"/"$log"
singularity exec "$kraken2" \
ktImportTaxonomy \
-q 2 \
-t 3 \
-tax "$CD"/databases/"$kronadb" \
"$TWD"/output/classified/${SLURM_ARRAY_TASK_ID}_kraken2_contigs_output.txt \
-o "$TWD"/output/classified/${SLURM_ARRAY_TASK_ID}_contigs_kronogram.html

# ITS identification and extraction with ITSx
mkdir -p "$TWD"/output/ITSx
name=$(echo $final_genome | cut -d'.' -f2 | rev | cut -d'/' -f1 | rev )
(date "+%T| Analysing ITS sequences using ITSx for $name") >> "$TWD"/"$log"
singularity exec "$itsx" \
ITSx \
-i "$final_genome" \
-o "$TWD"/output/ITSx/${SLURM_ARRAY_TASK_ID}_ITSx \
-t f \
--cpu "$cpus" \
# BLASTn to assign ITS with species name
singularity exec "$kraken2" \
blastn \
-db "$CD"/databases/fungi_ITS/ITS_RefSeq_Fungi \
-query "$TWD"/output/ITSx/${SLURM_ARRAY_TASK_ID}_ITSx.full.fasta \
-num_threads "$cpus" \
-outfmt "7 qacc ssciname sacc evalue" \
-out "$TWD"/output/ITSx/${SLURM_ARRAY_TASK_ID}_ITSx_nt.out
# extracts species name of top hit from BLASTn output
ITS_refseq=$(head -n6 output/ITSx/${SLURM_ARRAY_TASK_ID}_ITSx_nt.out | tail -nq | cut -f2 )
(date "+%T| Top hit for $name ITS: $ITS_refseq") >> "$TWD"/"$log"

# Quast evaluation
(date "+%T| Evaluating assembly ${SLURM_ARRAY_TASK_ID} with Quast") >> "$TWD"/"$log"
mkdir -p "$TWD"/output/quast/${SLURM_ARRAY_TASK_ID}
singularity exec "$qc" \
quast.py \
"$final_genome" \
-o "$TWD"/output/quast/${SLURM_ARRAY_TASK_ID} \
-t "$cpus" \
-L \
--fungus

quast_report="$TWD"/output/quast/${SLURM_ARRAY_TASK_ID}/report.tsv
if [ -f "$quast_report" ]; then
    (date "+%T| Quast evaluation of assembly ${SLURM_ARRAY_TASK_ID} completed") >> "$TWD"/"$log"
else
    (date "+%T| Quast evaluation of assembly ${SLURM_ARRAY_TASK_ID} did not complete correctly - check your inputs") >> "$TWD"/"$log"
    >&2 date "+%T| Quast evaluation of assembly ${SLURM_ARRAY_TASK_ID} did not complete correctly - check your inputs"
fi

# antiSMASH analysis
(date "+%T| Evaluating assembly ${SLURM_ARRAY_TASK_ID} with antiSMASH") >> "$TWD"/"$log"
mkdir -p "$TWD"/output/antismash/${SLURM_ARRAY_TASK_ID}
singularity exec "$antismash" \
run_antismash.py \
--cpus "$cpus" \
--taxon fungi \
--genefinding-tool glimmerhmm \
--fullhmmer \
--smcog-trees \
--cb-general \
--cb-subclusters \
--cb-knownclusters \
--minlength 5000 \
--output-dir "$TWD"/output/antismash/${SLURM_ARRAY_TASK_ID} \
"$final_genome"

antismash_report="$TWD"/output/antismash/${SLURM_ARRAY_TASK_ID}/${SLURM_ARRAY_TASK_ID}_scaffold.gbk
if [ -f "$antismash_report" ]; then
    (date "+%T| antiSMASH evaluation of assembly ${SLURM_ARRAY_TASK_ID} completed") >> "$TWD"/"$log"
else
    (date "+%T| antiSMASH evaluation of assembly ${SLURM_ARRAY_TASK_ID} did not complete correctly - check your inputs") >> "$TWD"/"$log"
    >&2 date "+%T| antiSMASH evaluation of assembly ${SLURM_ARRAY_TASK_ID} did not complete correctly - check your inputs" 
fi

(date "+%T| ~~ Assembly Evaluations Complete ~~") >> "$TWD"/"$log"
echo "" >> "$TWD"/"$log"

# removes unnecessary intermediate files

if [ "$space_saver" = "on" ]; then
    cd "$TWD"
    bytes=$(find -name "$SLURM_ARRAY_TASK_ID" -type d | du -sb | cut -f1)
    let size1=$bytes/1000000000
    (date "+%T| SpaceSaver is ON - Removing intermediate files") >> "$TWD"/"$log"
    rm -r "$TWD"/data/assembly/${SLURM_ARRAY_TASK_ID}/K*
    rm -r "$TWD"/data/assembly/${SLURM_ARRAY_TASK_ID}/tmp
    rm -r "$TWD"/data/assembly/${SLURM_ARRAY_TASK_ID}/misc
    rm -r "$TWD"/data/assembly/${SLURM_ARRAY_TASK_ID}/corrected
    rm "$TWD"/data/trimmed/"${SLURM_ARRAY_TASK_ID}_forward_trimmed_1U.fq.gz"
    rm "$TWD"/data/trimmed/"${SLURM_ARRAY_TASK_ID}_reverse_trimmed_2U.fq.gz" 
    rm "$TWD"/data/trimmed/"${SLURM_ARRAY_TASK_ID}_trimmomatic.log"
    bytes=$(find -name "$SLURM_ARRAY_TASK_ID" -type d | du -sb | cut -f1)
    let size2=$bytes/1000000000
    let size=size1-size2
    (date "+%T| SpaceSaver removed $size GB of unnecessary files") >> "$TWD"/"$log"
else
    (date "+%T| SpaceSaver is Off") >> "$TWD"/"$log"

fi

# computes the length of time script ran
end=$(date +%s)
runtime=$((end-start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))

echo "
Total runtime: $hours:$minutes:$seconds (hh:mm:ss)

If you found this script useful, or are having issues, send me an email at kelly.styles@vuw.ac.nz
Have a great day :D

-------------------------------- APPENDED SLURM LOG REPORT ---------------------------------------

" >> "$TWD"/"$log"

seff ${SLURM_JOB_ID}-${SLURM_ARRAY_TASK_ID} >> "$TWD"/"$log"
