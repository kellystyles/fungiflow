#!/usr/bin/env/python
import os
import sys
import re
import pandas as pd

"""
This script will find funannotate directories in the master directory and then
prepare a table containing boolean data about the presence/absence of the 
results file for each major step in the funannotate pipeline.

Usage `python3 collate_completion.py target_directory`
"""

print(f"Target directory : {sys.argv[1:]}")

target_dir = sys.argv[1:]
output_csv = os.path.join(target_dir, "completion_state.csv")

# Get the list of all directories in the folder
directories = [d for d in os.listdir(target_dir) if os.path.isdir(d)]

# Create an empty dataframe with columns for each file
columns = ['array', 'coverage', 'trimmed', 'kraken2', 'assembly', 
           'funannotate_predict', 'eggnog', 'funannotate_annotate', 
           'ITSx', 'antismash', 'quast'] 
df = pd.DataFrame(columns=columns) 
  
# Loop through all directories and check for the presence of each file 
for directory in directories: 
    # Boolean values to indicate presence/absence of files 
    coverage = os.stat(f'{directory}/trimmed/{directory}_trimmed_1P.fq.gz').st_size
    trimmed = os.path.exists(f'{directory}/trimmed/{directory}_trimmed_1P.fq.gz') and os.stat(f'{directory}/trimmed/{directory}_trimmed_1P.fq.gz').st_size > 0
    kraken2 = os.path.exists(f'{directory}/kraken2/{directory}_unclass_1.fq') and os.stat(f'{directory}/kraken2/{directory}_unclass_1.fq').st_size > 0
    assembly = os.path.exists(f'{directory}/assembly/scaffolds.fasta') and os.stat(f'{directory}/assembly/scaffolds.fasta').st_size > 0
    funannotate_predict = os.path.exists(f'{directory}/funannotate/predict_results/{directory}.gbk') and os.stat(f'{directory}/funannotate/predict_results/{directory}.gbk').st_size > 0
    eggnog = os.path.exists(f'{directory}/funannotate/eggnog/{directory}.emapper.annotations') and os.stat(f'{directory}/funannotate/eggnog/{directory}.emapper.annotations').st_size > 0
    funannotate_annotate = os.path.exists(f'{directory}/funannotate/annotate_results/{directory}.gbk') and os.stat(f'{directory}/funannotate/annotate_results/{directory}.gbk').st_size > 0
    ITSx = os.path.exists(f"{directory}/ITSx/{directory}.summary.txt") and os.stat(f"{directory}/ITSx/{directory}.summary.txt").st_size > 0
    antismash = os.path.exists(f'{directory}/antismash/{directory}.json') and os.stat(f'{directory}/antismash/{directory}.json').st_size > 0
    quast = os.path.exists(f'{directory}/quast/report.pdf') and os.stat(f'{directory}/quast/report.pdf').st_size > 0
    dict = {'array':directory, 'trimmed':trimmed, 'kraken2':kraken2, 'assembly':assembly, 
           'funannotate_predict':funannotate_predict, 'eggnog':eggnog, 'funannotate_annotate':funannotate_annotate, 
           'ITSx':ITSx, 'antismash':antismash, 'quast':quast}
    df = df.append(dict, ignore_index=True)

print(df)
df.to_csv(output_csv)
print(f"Results saved to {output_csv}")