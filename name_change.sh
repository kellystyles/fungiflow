#!/bin/bash

# Changes sequence reads filenames so they can be processed in a SLURM array
# Usage: bash name_change.sh directory

echo "Changing sequence filenames in the following directory: $1"

# loop to change forward/1 sequence reads
c=0
for i in $(ls -v $1/*1.f*.gz); do
    c=$(($c+1));
    A=$(echo $i | cut -d'/' -f3); # cut command only required if you want to remove a prefix; otherwise A=$(echo $i)
    echo " - $A renamed to $c" >> filenames.txt;
    mv $i "$c""_1.fq.gz";
done
# loop to change reverse/2 sequence reads
c=0
for i in $(ls -v $1/*2.f*.gz); do
    c=$(($c+1));
    A=$(echo $i | cut -d'/' -f3); # cut command only required if you want to remove a prefix; otherwise A=$(echo $i)
    echo " - $A renamed to $c" >> filenames.txt;
    mv $i "$c""_2.fq.gz";
done
