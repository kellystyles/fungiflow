#!/bin/bash

# Adds prefix to sequence read filenames so they can be processed in a SLURM array
# Usage: bash name_change.sh directory

echo "Changing sequence filenames in the following directory: $1"

dir=$(echo $1 | rev | cut -d"/" -f2- | rev)  # get directory path

for i in $(ls -v ${dir}/*_1*.f*); do
    c=$((${c}+1));      # iterator
    i1=$( basename $i )
    i2=$( basename $i1 1.fq.gz)2.fq.gz
    #echo $i1 $i2
    mv ${dir}/${i1} ${dir}/${c}_${i1};
    mv ${dir}/${i2} ${dir}/${c}_${i2};
    if [ -f ${dir}/${c}_${i1} ]; then
        echo " - ${dir}/${i1} renamed to ${dir}/${c}_${i1}"
        echo " - ${dir}/${i1} renamed to ${dir}/${c}_${i1}" >> filenames.txt;
    fi
    if [ -f ${dir}/${c}_${i2} ]; then
        echo " - ${dir}/${i2} renamed to ${dir}/${c}_${i2}"
        echo " - ${dir}/${i2} renamed to ${dir}/${c}_${i2}" >> filenames.txt;
    fi
done