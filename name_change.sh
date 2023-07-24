#!/bin/bash

# Adds prefix to sequence read filenames so they can be processed in a SLURM array
# Usage: bash name_change.sh directory

echo "Changing sequence filenames in the following directory: $1"

dir=$(echo $1 | rev | cut -d"/" -f2- | rev)  # get directory path

# loop to change sequence reads
c=0
for i in $(ls -v ${dir}/*[1F].f*.gz); do
    c=$((${c}+1));      # iterator
    F=$(echo ${i} | rev | cut -d'/' -f1 | rev); # get filename only
    echo " - ${i} renamed to ${dir}/${c}_${F}" >> filenames.txt;
    echo $i ${dir}/${c}_${F}
    mv ${i} ${dir}/${c}_${F};
    i2=${i//1/2} 
    i2=${i2//F/R} && echo $i2
    R=$(echo ${i2} | rev | cut -d'/' -f1 | rev)
    echo " - ${i2} renamed to ${dir}/${c}_${R}" >> filenames.txt;
    echo $i2 ${dir}/${c}_${R}
    mv ${i2} ${dir}/${c}_${R};
done

