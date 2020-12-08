#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --partition=quicktest
#SBATCH --time=01:00:00
#SBATCH --job-name=report_collator
#SBATCH -o /nfs/scratch/styleske/fungiflow/scripts/report_collator_%A_%a.out
#SBATCH -e /nfs/scratch/styleske/fungiflow/scripts/report_collator_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kelly.styles@vuw.ac.nz

#################################
### FUNGIFLOW REPORT COLLATOR ###
#################################

# This script collates all the fungiflow reports into a single meaningful summation report
#
# Usage:
# 1) edit the variable 'project' from this file to the name of your project
# 2) submit script to SLURM "sbatch report_collator.sh" from /yourpath/to/fungiflow/scripts/
#
# Outputs:
# - a table of all the Quast output in fungiflow/output/collated_quast_reports.tsv
# - a table counting the frequency of each biosynthetic gene cluster (BGC) from each antismash report;
#   in fungiflow/output/antismash/collated_antismash_results.tsv 

project="mixed"
collator="utility_v1.sif"

# makes logfile
CD=$(cd .. && pwd)
TWD="$CD"/projects/"$project" # working directory
time=$(date "+%H%M_%d_%m_%Y") # timestamp
log="$time"_"$project"_report_collator.txt
touch "$TWD"/"$log"

# Quast reports
echo "~~~ Report collation of $project FungiFlow output reports ~~~

The working directory is: $TWD

Date: $(date "+%A %d %b %Y")
Time: $(date "+%T")

####   ANALYSIS BEGINS    ####" >> "$TWD"/"$log"
(date "+%T| Collating Quast reports...") >> "$TWD"/"$log"

cd "$TWD"
awk -F '\t' '{_[FNR]=(_[FNR] OFS $2)}END{for (i=1; i<=FNR; i++) {sub(/^ /,"",_[i]); print _[i]}}' \
output/quast/*/report.tsv > tmp1 && tr ' ' '\t' < tmp1 > tmp2
find -name "report.tsv" > quastreportlist
x=$(head -1 quastreportlist)
cut -f1 $x > tmp1
paste tmp1 tmp2 > output/quast/collated_quast_reports.tsv
rm tmp* quastreportlist

(date "+%T| Finished collating Quast reports") >> "$TWD"/"$log"

# antiSMASH reports
(date "+%T| Collating antiSMASH reports...") >> "$TWD"/"$log"

cd "$TWD"
find -name "*.gbk" | egrep "*scaffold*" > as_tmp1

module load singularity/3.6.1
singularity exec ../../images/"$collator" \
python3 ../../scripts/antismash_gbk_parser.py

mv "$TWD"/collated_antismash_results.tsv "$TWD"/output/antismash/collated_antismash_results.tsv
rm as_tmp1

(date "+%T| Finished collating antiSMASH reports...") >> "$TWD"/"$log"

echo "
If you found this script useful, or are having issues, send me an email at kelly.styles@vuw.ac.nz
Have a great day :D"

