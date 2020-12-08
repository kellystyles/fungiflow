#!/bin/bash
#SBATCH --cpus-per-task=24
#SBATCH --mem=125G
#SBATCH --partition=bigmem
#SBATCH --time=4-00:00:00
#SBATCH --job-name=itsx_v1
#SBATCH -o /nfs/scratch/styleske/fungiflow/scripts/itsx_%A_%a.out
#SBATCH -e /nfs/scratch/styleske/fungiflow/scripts/itsx_%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kelly.styles@vuw.ac.nz

project="mixed"

# makes logfile
CD=$(cd .. && pwd)
TWD="$CD"/projects/"$project" # working directory
time=$(date "+%H%M_%d_%m_%Y") # timestamp
log="$time"_"$project"_ITSx.txt
touch "$TWD"/"$log"

itsx="$CD/images/itsx_v1.sif"
kraken2_v1="$CD/images/kraken2_v1.sif"
export BLASTDB=$BLASTDB:$CD/databases/fungi_ITS/

cd "$TWD"
module load singularity/3.6.1
mkdir -p "$TWD"/output/ITSx

echo "~~~ ITSx analysis of $project assemblies ~~~"

find . -name "*_scaffold.fasta" > assemblies

while read -r line; do
    name=$( echo $line | cut -d'.' -f2 | rev | cut -d'/' -f1 | rev )
    if [ ! -f output/ITSx/"$name"_ITSx.full.fasta ]; then
        echo "Extracting ITS region for $name" >> "$TWD"/"$log"
        singularity exec "$itsx" \
        ITSx \
        -i "$line" \
        -o output/ITSx/"$name"_ITSx \
        -t f \
        --cpu 12
    fi
    if [ ! -f output/ITSx/"$name"_ITSx_nt.out ]; then    
        echo "Identifying ITS using BLASTn for $name" >> "$TWD"/"$log"
        singularity exec "$kraken2_v1" \
        blastn \
        -db "$CD"/databases/fungi_ITS/ITS_RefSeq_Fungi \
        -query output/ITSx/"$name"_ITSx.full.fasta \
        -num_threads 24 \
        -outfmt "7 ssciname sacc evalue" \
        > output/ITSx/"$name"_ITSx_nt.out
    fi
done < assemblies

echo "Parsing BLASTn reports..." >> "$TWD"/"$log"
find -name "*_ITSx_nt.out" > itsx_blastn_results.txt
while read -r line; do
    name=$(echo $line | rev | cut -d'/' -f1 | rev )
    ITS_refseq=$( head -n10 $line | tail -n5 )
    echo "$name
    $ITS_refseq" >> "$project"_its_results.tsv
done < itsx_blastn_results.txt
