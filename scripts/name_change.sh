#!/bin/bash

# Change filenames so they can be processed in a SLURM array
# Change the directory and file names for your files

c=0
for i in $(ls -v *_R1.f*q.gz); do
    c=$(($c+1));
    A=$(echo $i | cut -d'/' -f3); # cut command only required if you want to remove a prefix; otherwise A=$(echo $i)
    echo "$c""+""$A" >> filenames.txt;
    mv $i "$c""_R1.fq.gz";
done
c=0
for i in $(ls -v *R2.f*q.gz); do
    c=$(($c+1));
    A=$(echo $i | cut -d'/' -f3); # cut command only required if you want to remove a prefix; otherwise A=$(echo $i)
    echo "$c""+""$A" >> filenames.txt;
    mv $i "$c""_R2.fq.gz";
done
