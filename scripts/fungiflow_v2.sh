#!/bin/bash
#SBATCH --array=29-142%4
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=3G
#SBATCH --constraint=Intel
#SBATCH --partition=bigmem
#SBATCH --time=1-00:00:00
#SBATCH --job-name=fungiflow_v2
#SBATCH -o /nfs/scratch/styleske/fungiflow/scripts/fungiflow_v2_%A_%a.out
#SBATCH -e /nfs/scratch/styleske/fungiflow/scripts/fungiflow_v2_%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kelly.styles@vuw.ac.nz

# ----------------------------------------
# Welcome to the Fungiflow Analysis script
# ----------------------------------------
#
# This script assembles and evaluates short (Illumina) reads. 
# Firstly FastQC is used to evaluate reads, before parsing the reads for adapter sequences, 
# then trimming the reads using Trimmomatic. The trimmed reads are filtered using Kraken2 and
# the unclassified sequences (those that do not match anything in the standard database) preserved 
# and re-evaluated with FastQC. MultiQC will collate all the FastQC reports for easy viewing. 
# The cleaned reads are assembled using SPADES (or metaSPADES for metagenomic samples).
# Should this assembly fail, MaSurCa will be used to try to assemble the reads instead.
# The resulting scaffolds will be analysed with Quast and antiSMASH.
# Kraken2 will again be used to identify fungal sequences in assembly and the results
# will be plotted with Krona.
#
# The script will exit prematurely if the trimming, initial classification, or assembly steps fail.
#
# The FastQC emulator 'Falco' is also available for use in place of FastQC. Falco is a bit
# faster than FastQC so may be the better choice for larger datasets.
#
# Inputs:
# - raw paired Illumina reads in 'fq.gz' format (no pre-processed samples please)
# - make sure all the variables are set correctly
# - make sure that the QC, Kraken2, Assembly, and antiSMASH Singularity images are
#   in the images directory
# - make sure the kraken2 standard, kraken2 fungi, and krona taxonomy databses are
#   in the databases directory
#
# Please ensure that all raw reads fed into this script are prefixed with a number 
# and underscore like so: "X_" Otherwise, the SLURM task manager won't recognise your files
# To easily do this to a batch of files,'bash' the name_change.sh script in the 'scripts' folder.
# In addition, make sure that the array value in the SLURM header is correct and make sure that
# the cpu and memory values in the SLURM header match those in the 'Variables' section of this script.
# 
# Output:
# - the log file for each array task can be found in /project_directory
# - adapters can be found in /project_directory/data/adapters/ in 'fasta' format
# - trimmed reads can be found in /project_directory/data/trimmed/ in 'fg.gz' format
# - the FastQC and MultiQC output files can be found in /project_directory/output/(raw_QC and trimmed_QC)/
# - Kraken2-filtered reads can be found in /project_directory/data/classified/
#       - the classified reads are named <SLURM_ARRAY_TASK_ID>_cseqs_1.fq and <SLURM_ARRAY_TASK_ID>_cseqs_2.fq
#       - the unclassified reads are named <SLURM_ARRAY_TASK_ID>_useqs_1.fq and <SLURM_ARRAY_TASK_ID>_useqs_2.fq
# - the assembled reads can be found in /project_directory/data/assembly/<SLURM_ARRAY_TASK_ID>/
# - Quast output can be found in /project_directory/output/quast/${SLURM_ARRAY_TASK_ID}/
# - antiSMASH output can be found in /project_directory/output/antismash/${SLURM_ARRAY_TASK_ID}/
#
# Testing on a set of paired Illumina MiSeq reads of ~220 Mb each took about 2 hours with 24 cpus and 25 GB memory,
# although this will change depending on the amount of fungal reads recovered after filtering with Kraken2. 
# Feel free to play around with the computational specifications to find suitable values for your dataset. 
# I would recommend setting the initial job with a single array value to test the entire pipeline is working.
# If space is an issue, turn 'on' the space saver variable, which will delete various intermediate files.


### Variables ###
project="endolichenic_isolates"         # loads project name to target your project directory
type="isolate"      # what type of sample data are you assembling ('isolate' or a 'metagenomic')
cpus="24"               # Number of cpus to use per task
mem="72"                # Amount of memory per task in GB
eval="fastqc"           # read evaluator ('fastqc' or 'falco')
paired="yes"            # are the reads paired?
stddb="kraken2_standard" # the Kraken2 standard database
fundb="kraken2_fungi"   # the Kraken2 fungi database
kronadb="krona"         # the Krona taxonomy database name
sif1="qc_v1.sif"        # Singularity image to use for QC
sif2="kraken2_v1.sif"   # Singularity image to use for Kraken2
sif3="assembly_v2.sif"  # Singularity image to use for Assembly
sif4="antismash_5_1_ver3" # Singularity image to use of antiSMASH
space_saver="on"        # removes intermediate files ('on' or 'off')

###############################
####     Script Begins     ####
###############################

start=$(date +%s)

# makes logfile
CD=$(cd .. && pwd)
TWD="$CD"/projects/"$project" # working directory
time=$(date "+%H%M_%d_%m_%Y") # timestamp
log="$time"_"$project"_fungiflow_${SLURM_ARRAY_TASK_ID}.txt
touch "$TWD"/"$log"

# Initialize Singularity module on Raapoi
module load singularity/3.6.1
qc="$CD"/images/"$sif1"
kraken2="$CD"/images/"$sif2"
assembly="$CD"/images/"$sif3"
antismash="$CD"/images/"$sif4"

echo "~~~ Analysis of $project ${SLURM_ARRAY_TASK_ID} Illumina reads ~~~

The working directory is: $TWD

Job name:    $SLURM_JOB_NAME - $project
Job ID:      $SLURM_JOB_ID
Job array:   $SLURM_ARRAY_TASK_ID 

Images:
 - QC           = $qc
 - Kraken2      = $kraken2
 - Assembly     = $assembly
 - antiSMASH    = $antismash 

Date: $(date "+%A %d %b %Y")
Time: $(date "+%T")
Software versions:
 - FastQC v0.11.9
 - Falco v0.2.1
 - Trimmomatic v0.39
 - MultiQC v1.9
 - Kraken2 v2.1.0
 - SPADES v3.14.1
 - Krona v2.7.1

Files to be processed:" >> "$TWD"/"$log"
file_list=$(ls -v "$TWD"/data/raw/*.fq*)
for line in $file_list; do
    file=${line##*/}
    [[ "$file" = ${SLURM_ARRAY_TASK_ID}_* ]] && \
    echo " - $file" >> "$TWD"/"$log"
done

echo "
####   ANALYSIS BEGINS    ####

       Quality Control
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" >> "$TWD"/"$log"

# ```
# QC
# ```

# fastqc on raw reads
mkdir -p "$TWD"/output/raw/
(date "+%T| Performing fastqc on ${SLURM_ARRAY_TASK_ID} raw paired reads") >> "$TWD"/"$log"
singularity exec "$qc" \
"$eval" \
-o "$TWD"/output/raw/ \
"$TWD"/data/raw/${SLURM_ARRAY_TASK_ID}_*.fq*

# uses adap_ID.sh to find all adapters in each file and writes them to a file
(date "+%T| Identifying adapter sequences in ${SLURM_ARRAY_TASK_ID} reads") >> "$TWD"/"$log"
mkdir -p "$TWD"/data/adapters/
singularity exec "$qc" \
bash /adap_ID.sh \
"$TWD"/data/raw/${SLURM_ARRAY_TASK_ID}_R1.fq.gz | grep "^>" -A 1 --no-group-separator > \
"$TWD"/data/adapters/${SLURM_ARRAY_TASK_ID}_adapter_file.fasta
singularity exec "$qc" \
bash /adap_ID.sh \
"$TWD"/data/raw/${SLURM_ARRAY_TASK_ID}_R2.fq.gz | grep "^>" -A 1 --no-group-separator >> \
"$TWD"/data/adapters/${SLURM_ARRAY_TASK_ID}_adapter_file.fasta
    
# appends TruSeq2-PE adapters onto end of created adapter files
singularity exec "$qc" \
cat /TruSeq3-PE-2.fa >> "$TWD"/data/adapters/${SLURM_ARRAY_TASK_ID}_adapter_file.fasta

# trims the reads with Trimmomatic
(date "+%T| Trimming Illumina ${SLURM_ARRAY_TASK_ID} reads with Trimmomatic") >> "$TWD"/"$log"
mkdir -p "$TWD"/data/trimmed/
singularity exec "$qc" \
TrimmomaticPE \
-threads "$cpus" \
-trimlog "$TWD"/data/trimmed/${SLURM_ARRAY_TASK_ID}_trimmomatic.log \
"$TWD"/data/raw/${SLURM_ARRAY_TASK_ID}_R1.fq* "$TWD"/data/raw/${SLURM_ARRAY_TASK_ID}_R2.fq* \
"$TWD"/data/trimmed/${SLURM_ARRAY_TASK_ID}_forward_trimmed_1P.fq.gz \
"$TWD"/data/trimmed/${SLURM_ARRAY_TASK_ID}_forward_trimmed_1U.fq.gz \
"$TWD"/data/trimmed/${SLURM_ARRAY_TASK_ID}_reverse_trimmed_2P.fq.gz \
"$TWD"/data/trimmed/${SLURM_ARRAY_TASK_ID}_reverse_trimmed_2U.fq.gz \
ILLUMINACLIP:"$TWD"/data/adapters/${SLURM_ARRAY_TASK_ID}_adapter_file.fasta:2:30:10:4:4:/true \
TRAILING:9 SLIDINGWINDOW:4:15 MINLEN:36

trimmedf="$TWD"/data/trimmed/${SLURM_ARRAY_TASK_ID}_forward_trimmed_1P.fq.gz
trimmedr="$TWD"/data/trimmed/${SLURM_ARRAY_TASK_ID}_reverse_trimmed_2P.fq.gz
if [ -f "$trimmedf" ] && [ -f "$trimmedr" ]; then
    (date "+%T| Trimming ${SLURM_ARRAY_TASK_ID} reads with Trimmomatic completed") >> "$TWD"/"$log"
else
    (date "+%T| Trimming ${SLURM_ARRAY_TASK_ID} reads with Trimmomatic did not complete correctly - check your inputs") >> "$TWD"/"$log"
    >&2 date "+%T| Trimming ${SLURM_ARRAY_TASK_ID} reads with Trimmomatic did not complete correctly - check your inputs"
    exit 1
fi

# Kraken2 filtering of trimmed reads
(date "+%T| Beggining Kraken2 analysis of file for array ${SLURM_ARRAY_TASK_ID}") >> "$TWD"/"$log"
mkdir -p "$TWD"/data/classified/
mkdir -p "$TWD"/output/classified/
singularity exec "$kraken2" \
kraken2 \
--db $CD/databases/$stddb \
--threads "$cpus" \
--use-names \
--paired \
--output "$TWD"/output/classified/${SLURM_ARRAY_TASK_ID}_kraken2_reads_output.txt \
--classified-out "$TWD"/data/classified/${SLURM_ARRAY_TASK_ID}_cseqs#.fq \
--unclassified-out "$TWD"/data/classified/${SLURM_ARRAY_TASK_ID}_useqs#.fq \
"$TWD"/data/trimmed/${SLURM_ARRAY_TASK_ID}_forward_trimmed_1P.fq.gz \
"$TWD"/data/trimmed/${SLURM_ARRAY_TASK_ID}_reverse_trimmed_2P.fq.gz

kraken2f="$TWD"/data/classified/${SLURM_ARRAY_TASK_ID}_useqs_1.fq
kraken2r="$TWD"/data/classified/${SLURM_ARRAY_TASK_ID}_useqs_2.fq
if [ -f "$kraken2f" ] && [ -f "$kraken2r" ]; then
    (date "+%T| Kraken2 analysis of ${SLURM_ARRAY_TASK_ID} trimmed reads completed") >> "$TWD"/"$log"
else
    (date "+%T| Kraken2 analysis of ${SLURM_ARRAY_TASK_ID} trimmed reads did not complete correctly - check your inputs") >> "$TWD"/"$log"
    >&2 date "+%T| Kraken2 analysis of ${SLURM_ARRAY_TASK_ID} trimmed reads did not complete correctly - check your inputs"
    exit 1
fi

# fastqc on trimmed and filtered reads
mkdir -p "$TWD"/output/trimmed/
(date "+%T| Performing fastqc on ${SLURM_ARRAY_TASK_ID} trimmed paired reads") >> "$TWD"/"$log"
singularity exec "$qc" \
"$eval" \
-o "$TWD"/output/trimmed/ \
"$TWD"/data/classified/${SLURM_ARRAY_TASK_ID}_useqs_*.fq

(date "+%T| ~~ QC Complete ~~") >> "$TWD"/"$log"

# ```
# Assembly
# ```

echo "
           Assembly 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" >> "$TWD"/"$log"

if [ "$type" = "metagenomic" ]; then
    assembler="MetaSPADES"
    # Assembly with MetaSPADES
    (date "+%T| Assembling trimmed ${SLURM_ARRAY_TASK_ID} reads with MetaSPADES ") >> "$TWD"/"$log"
    mkdir -p "$TWD"/data/assembly/${SLURM_ARRAY_TASK_ID}
    singularity exec "$assembly" \
    metaspades.py \
    --threads "$cpus" \
    --memory "$mem" \
    -k 21,33,55,77,99,127 \
    --pe1-1 "$TWD"/data/classified/${SLURM_ARRAY_TASK_ID}_useqs_1.fq \
    --pe1-2 "$TWD"/data/classified/${SLURM_ARRAY_TASK_ID}_useqs_2.fq \
    -o "$TWD"/data/assembly/${SLURM_ARRAY_TASK_ID}
    (date "+%T| Finished assembling ${SLURM_ARRAY_TASK_ID} reads with MetaSPADES") >> "$TWD"/"$log"
else
    assembler="SPADES"
    # Assembly with SPADES
    (date "+%T| Assembling trimmed ${SLURM_ARRAY_TASK_ID} reads with SPADES ") >> "$TWD"/"$log"
    mkdir -p "$TWD"/data/assembly/${SLURM_ARRAY_TASK_ID}
    singularity exec "$assembly" \
    spades.py \
    --threads "$cpus" \
    --memory "$mem" \
    --careful \
    -k 21,33,55,77,99,127 \
    --pe1-1 "$TWD"/data/classified/${SLURM_ARRAY_TASK_ID}_useqs_1.fq \
    --pe1-2 "$TWD"/data/classified/${SLURM_ARRAY_TASK_ID}_useqs_2.fq \
    -o "$TWD"/data/assembly/${SLURM_ARRAY_TASK_ID}
    (date "+%T| Finished assembling ${SLURM_ARRAY_TASK_ID} reads with SPADES") >> "$TWD"/"$log"
fi

mv "$TWD"/data/assembly/${SLURM_ARRAY_TASK_ID}/scaffolds.fasta "$TWD"/data/assembly/${SLURM_ARRAY_TASK_ID}/${SLURM_ARRAY_TASK_ID}_scaffold.fasta

final_genome="$TWD"/data/assembly/${SLURM_ARRAY_TASK_ID}/${SLURM_ARRAY_TASK_ID}_scaffold.fasta
if [ -f "$final_genome" ]; then
    (date "+%T| "$assembler" assembly of ${SLURM_ARRAY_TASK_ID} trimmed reads completed") >> "$TWD"/"$log"
else
    (date "+%T| "$assembler" assembly of ${SLURM_ARRAY_TASK_ID} trimmed reads did not complete correctly") >> "$TWD"/"$log"
    >&2 date "+%T| "$assembler" assembly of ${SLURM_ARRAY_TASK_ID} trimmed reads did not complete correctly"
    
    (date "+%T| Trying MaSurCa assembly of ${SLURM_ARRAY_TASK_ID} trimmed reads") >> "$TWD"/"$log"
    mkdir -p "$TWD"/data/assembly/${SLURM_ARRAY_TASK_ID}/masurca/

#### generating masurca config file ####
echo "DATA
PE= BU 250 50 $TWD/$directory/${SLURM_ARRAY_TASK_ID}_useqs_1.fq $TWD/$directory/${SLURM_ARRAY_TASK_ID}_useqs_2.fq
END

PARAMETERS
EXTEND_JUMP_READS=0
GRAPH_KMER_SIZE = auto
USE_LINKING_MATES = 1
LHE_COVERAGE=30
MEGA_READS_ONE_PASS=0
LIMIT_JUMP_COVERAGE = 300
CA_PARAMETERS =  cgwErrorRate=0.15
NUM_THREADS = $cpus
JF_SIZE = 800000000
SOAP_ASSEMBLY=0
END" > "$TWD"/data/assembly/${SLURM_ARRAY_TASK_ID}/masurca/${SLURM_ARRAY_TASK_ID}_masurca_config

#### end of config file ####

    singularity exec "$assembly" \
    masurca "$TWD"/data/assembly/${SLURM_ARRAY_TASK_ID}/masurca/${SLURM_ARRAY_TASK_ID}_masurca_config
    singularity exec "$assembly" \
    bash "$TWD"/data/assembly/${SLURM_ARRAY_TASK_ID}/masurca/assemble.sh
    
    mv "$TWD"/data/assembly/${SLURM_ARRAY_TASK_ID}/masurca/CA/final.genome.scf.fasta "$TWD"/data/assembly/${SLURM_ARRAY_TASK_ID}/masurca/CA/${SLURM_ARRAY_TASK_ID}_scaffold.fasta
    final_genome="$TWD"/data/assembly/${SLURM_ARRAY_TASK_ID}/masurca/CA/${SLURM_ARRAY_TASK_ID}_scaffold.fasta
    if [ -f "$final_genome" ]; then
        (date "+%T| MaSurCa assembly of ${SLURM_ARRAY_TASK_ID} trimmed reads complete") >> "$TWD"/"$log"
    else
        (date "+%T| MaSurCa assembly of ${SLURM_ARRAY_TASK_ID} trimmed reads did not complete correctly - check your inputs") >> "$TWD"/"$log"
        >&2 date "+%T| MaSurCa assembly of ${SLURM_ARRAY_TASK_ID} trimmed reads did not complete correctly - check your inputs"
        exit 1
    fi
fi

(date "+%T| ~~ Assembly Complete ~~") >> "$TWD"/"$log"

# ```
# Assembly Evaluations
# ```

echo "
     Assembly Evaluations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" >> "$TWD"/"$log"

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
"$TWD"/output/classified/${SLURM_ARRAY_TASK_ID}_kraken2_contigs_report.txt \
-o "$TWD"/output/classified/${SLURM_ARRAY_TASK_ID}_contigs_kronogram.html

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

vuw-job-report "$SLURM_JOB_ID"_${SLURM_ARRAY_TASK_ID} >> "$TWD"/"$log"
