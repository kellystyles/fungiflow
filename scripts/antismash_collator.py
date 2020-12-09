#!/usr/bin/python

import os
import os.path
import json
import subprocess
import argparse

version = 0.1

parser = argparse.ArgumentParser(description="""\
----------------------- antiSMASH output collator v{version} ----------------------- \
\
This script will find all antiSMASH '.json' output files from a given 'input_directory' (including subdirectories), \
parses and then outputs a combined '.json' file with the 'output_name' as output. \
\
~~ JSON parse functions written by Jeremy Owen | \
antismash_collator.py is maintained by Kelly Styles - kelly.styles@vuw.ac.nz ~~ \
""".format(version=version))

# positional arguments                                 
parser.add_argument("input_directory", 
                    help = "path to input directory")
parser.add_argument("output_name", 
                    help="name for results")

# optional arguments
parser.add_argument("-c", "--cpus",  
                    default=2, 
                    type=int,
                    help="number of cpus to use, default=2")
parser.add_argument("-t", "--types",
                    type=set,
                    help="target cluster types as a python set")
parser.add_argument("-p", "--partial",  
                    default="False",
                    action='store_true', 
                    help="partial clusters (on contig edge) will be included")
parser.add_argument("-taxa", "--taxonomy",  
                    default="bacteria",
                    help="designate which BGCs to search for based on taxonomy ('bacteria' or 'fungi')")
parser.add_argument("-r", "--rebuild",  
                    default="False",
                    help="rebuilds collated JSON file with designated antiSMASH Singularity image (/path/to/image)")
args = parser.parse_args()

def ParseOne(file, 
             cluster_types, 
             partial=False):
    
    """
    Parses a single JSON file and identifies cluster types found in the list (default = all cluster types)
    If allow partial is True, partial clusters (on contig edge) will be included. output_name is the name of the folder
    the json file is found in
    Input: '.json' file
    Output: Returns a tuple of two lists:
        tup[0] = list of record (containing features etx)
        tup[1] = corresponding 'timings' entries
    """
    ID=0
    with open(file) as F:
        parsed = json.load(F)
    
    output_name = file.split("/")[-1] # takes name from folder the .json file is in
    timings_to_return = {}
    records_to_return = []
    
    for i in range(len(parsed['records'])):
        for j in range(len(parsed['records'][i]['features'])):
            record = parsed['records'][i]['features'][j]
            
            if ('cand_cluster' in record['type']) \
            and (partial or 
                 'False' in record["qualifiers"]["contig_edge"]):
                
                is_target = False
                for index in range(len(record["qualifiers"]['product'])):
                    if record["qualifiers"]['product'][index] in cluster_types:
                        is_target = True
                        
                if is_target:
                    cluster_ID = output_name + "_" + str(ID)
                    record_to_add = parsed['records'][i]
                    timing_to_add = parsed['timings'][record_to_add['id']]
                    record_string = json.dumps(record_to_add).replace(record_to_add['id'],cluster_ID)
                    record_to_add = json.loads(record_string)
                    timing_string = json.dumps(timing_to_add).replace(record_to_add['id'],cluster_ID)
                    timing_to_add = json.loads(timing_string)
                                                        
                    records_to_return.append(record_to_add)
                    timings_to_return[cluster_ID] = timing_to_add
                    ID += 1
    
    return (records_to_return,timings_to_return)

def GetJsonPaths(target_dir):
    """
    Finds all '.json' files in a directory and its subdirectories
    Returns a list of filepaths for '.json' files
    """
    json_paths = []
    for dirpath, dirnames, filenames in os.walk(target_dir):
        for filename in [f for f in filenames if f.endswith(".json")]:
            json = os.path.join(dirpath, filename)
            json_paths.append(json)
    return json_paths

def ParseFolder(target_dir, 
                types, 
                partial = False, 
                antismash_version = '5.1.0', 
                input_file= "multi", 
                taxon = 'bacteria', 
                schema = 1):
    """
    Parses a folder containing antismash outputs (subdirectories)
    Runs ParseOne on each antismash output
    Returns a '.json' file containing all BGCs of the designated cluster types
    """
    new_json = {'version':antismash_version, 
                'input_file':input_file, 
                'records':[], 
                'timings':{}, 
                'taxon':taxon, 
                'schema': schema}
    json_paths = GetJsonPaths(target_dir)
    total = len(json_paths)
    current = 1
    for json in json_paths:
        print("processing file {} of {}".format(current, total))
        print("file name: " + json)
        tup = ParseOne(json,cluster_types=types,partial=partial)
        new_json['records'].extend(tup[0])
        for key in tup[1]:
            new_json['timings'][key] = tup[1][key]
        
        current += 1
    
    return new_json

def jsonDump(parsed_folder, output_name):
    """
    Takes the output of ParseFolder function and writes it to a new json file
    Input: output of ParseFolder function
    Output: ".json" folder
    """
    with open(output_name, 'w') as F:
        json.dump(parsed_folder,F)
        
def generate_bgc_html(image_path, json, output_name, cpus):
    """
    Rebuilds the antismash index.html file from a '.json' file
    input: '.json' file, relative path of antismash Singularity image
    output: rebuilt antismash directory
    """
    cmd1 = "sbatch --cpus-per-task=" + str(cpus) + " antismash_rebuild.sh " + image_path + " " + json + " " + output_name
    
    print("Rebuilding multijson file with antiSMASH...")
    subprocess.call(cmd1, shell=True)

__version__ = "0.1"

# Begin actual script

print("----------------------- antiSMASH report collator v{version} -----------------------".format(version=version))
target_dir = args.input_directory
print("Input folder is at %s" % target_dir)
output_name = args.output_name
print("Output name is %s" % output_name)

if args.cpus:
    cpus = args.cpus
    print("%s cpus will be used" % cpus)
if args.types:
    types = args.types
    print("The following cluster types will be filtered: %s" % types)
else:
    types = {"T1PKS", "T2PKS", "T3PKS", "transAT-PKS", "transAT-PKS-like", "PpyS-KS", \
            "hglE-KS", "CDPS", "PKS-like", "arylpolyene", "resorcinol", "ladderane", \
            "PUFA", "nrps", "nrps-like", "thioamide-NRP", "terpene", "lanthipeptide", \
            "lipolanthine", "bacteriocin", "betalactone", "thiopeptide", "linaridin", \
            "cyanobactin", "glycocin", "LAP", "lassopeptide", "sactipeptide", "bottromycin", \
            "head_to_tail", "microviridin", "proteusin", "blactam", "amglyccycl", "aminocoumarin", \
            "siderophore", "ectoine", "butyrolactone", "indole", "nucleoside", "phosphoglycolipid", \
            "melanin", "oligosaccharide", "furan", "hserlactone", "phenazine", "phosphonate", \
            "fused", "PBDE", "acyl_amino_acids", "tropodithietic-acid", "NAGGN", "RaS-RiPP", \
            "fungal-RiPP", "TfuA-related", "other", "saccharide", "fatty_acid", "halogenated"}
if args.partial:
    partial = args.partial
    print("Allow partial clusters is set to %s" % partial)
if args.taxonomy:
    taxonomy = args.taxonomy
    print("The taxonomy is %s" % taxonomy)
if args.rebuild is not None:
    antismash_image = args.rebuild
    print("The collated JSON file will be rebuilt with antiSMASH using %s" % antismash_image)

output_file = output_name + ".json"
print("New JSON file '%s' has been made" % output_file)
parsed_folder = ParseFolder(target_dir, types, partial, taxonomy)
bgc_count = len(parsed_folder['records'])
print("Found %s BGCs" % bgc_count)
print("Writing BGCs to new JSON file...")
jsonDump(parsed_folder, output_file)
if args.rebuild is not None:
    generate_bgc_html(antismash_image, output_file, output_name, cpus)
