import os
import argparse
import subprocess
import sys
import pyhmmer
import multiprocessing
import ast
import pandas as pd
from Bio import SeqIO

version = 1
f"""
                    |-|-|-|- Cluster Search v1 -|-|-|-|

This script will pull out sequence data containing hits to input pHMM models 
from a GBK file. Essentially, you can input pHMMs for your gene of interest and
then it will collate all top hits for each pHMM in a new 'cluster' GBK file. 

Usage: python3 cluster_search_v{str(version)}.py -i gbk_dir -d phmm_dir -r required_models

Input:
    --gbk_dir   =   directory containing genbank files
    --phmm_dir  =   directory containing pHMM files
    --required  =   (optional) pHMMs that must be included in output clusters

                    |-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|
"""

def get_args():
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser(
        description="Finds and extracts clusters of biosynthetic genes in an annotated genome using given HMM models")
    parser.add_argument('-g', '--gbk_dir', action='store',
                        help='genbank file(s)', required=True)
    parser.add_argument('-p', '--phmm_dir', action='store',
                        help='directory with pHMMs', required=True)
    parser.add_argument('-a', '--alignment', action='store',
                        help='Multiple sequence alignment to build pHMM from', required=False)
    parser.add_argument('-r', '--required', action='store',
                        help='pHMMs required to be in the cluster')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        exit(1)
    
    return parser.parse_args()

def check_versions():
    """
    Print installed software versions
    """
    #print("Installed Software:\n")

    # Python modules
    #print(f"Python: {sys.version}")
    #print(f"BioPython: {Bio.__version__}")
    #print(f"pandas: {pd.__version__}")
    
    # Bash modules
    cmds = ["blastp -h", "hmmsearch -h", "muscle"]
    for cmd in cmds:
        #print(cmd)
        subprocess.run([cmd], shell=True)

def get_cds(gbk, out_filename):
    """
    Extract protein sequences from a GBK file and write to a multi-FASTA file
    """
    if os.path.exists("cds") is False:
        os.makedirs("cds", exist_ok=True)
    if not os.path.isfile(out_filename):
        with open(gbk) as handle, open(out_filename, "w") as output_handle:
            for seq_record in SeqIO.parse(handle, "genbank"):
                for seq_feature in seq_record.features:
                    if seq_feature.type == "CDS":
                        try:
                            assert len(seq_feature.qualifiers['translation']) == 1
                            output_handle.write(">%s \n%s\n" % (
                                seq_feature.qualifiers['protein_id'][0],     # change protein_id to appropriate name value, depending on GBK source
                                seq_feature.qualifiers['translation'][0]))
                        except KeyError:
                            continue
    return out_filename

def build_hmm(msa, output_file):
    """
    Builds a binary HMM 'h3m' file using pyhmmer from an alignment file
    """
    alphabet = pyhmmer.easel.Alphabet.amino()
    builder = pyhmmer.plan7.Builder(alphabet)
    background = pyhmmer.plan7.Background(alphabet)
    
    with pyhmmer.easel.MSAFile(msa, digital=True, alphabet=alphabet) as msa_file:
        msa_d = msa_file.read()
    prefix = msa.split(".")[0]
    msa_d.name = prefix.encode()
    hmm, _, _ = builder.build_msa(msa_d, background)

    with open(output_file, "wb") as f:
        print(f"Writing HMM to {output_file}")
        hmm.write(f)
    
    return msa_d

def read_hmm(h3m_file):
    with pyhmmer.HMMFile(h3m_file) as hmm_file:
        hmms = list(hmm_file)
    
def multi_hmmsearch(hmm_list, fasta):
    from typing import Iterable
    from pyhmmer.plan7 import HMMFile
    from pyhmmer.easel import Bitfield, TextSequence
    from pyhmmer.hmmer import hmmsearch
    import collections

    alphabet = pyhmmer.easel.Alphabet.amino()

    with pyhmmer.easel.SequenceFile(fasta, digital=True) as seqs_file:
        proteins = seqs_file.read_block()
    print(f"Loaded {len(proteins)} protein CDS from {fasta}")
        
    loaded = []
    for f in hmm_list:
        with pyhmmer.plan7.HMMFile(f) as hmm_file:
            hmm = hmm_file.read()
            loaded.append(hmm)
        
    Result = collections.namedtuple("Result", ["query", "hmm", "bitscore"])
    results = []
    for hits in pyhmmer.hmmsearch(loaded, proteins):
        hmm = hits.query_name.decode()
        for hit in hits:
            if hit.included:
                results.append(Result(hit.name.decode(), hmm, hit.score))
                #print(results)
    
    return results

def trusted_hmmsearch(hmms, fasta, trusted_cutoff_file):
    from typing import Iterable
    from pyhmmer.plan7 import HMMFile
    from pyhmmer.easel import Bitfield, TextSequence
    from pyhmmer.hmmer import hmmsearch
    import collections

    alphabet = pyhmmer.easel.Alphabet.amino()

    cutoffs = {}
    with open(trusted_cutoff_file, "r") as tc_file:
        for line in tc_file:
            line = line.split()
            cutoffs[line[0]] = float(line[1])
    
    with pyhmmer.easel.SequenceFile(fasta, digital=True) as seqs_file:
        proteins = seqs_file.read_block()
    print(f"Loaded {len(proteins)} protein CDS from {fasta}")
    
    loaded = []
    for f in hmms:
        with pyhmmer.plan7.HMMFile(f) as hmm_file:
            hmm = hmm_file.read()   
        cutoff = cutoffs[hmm.name.decode()]    
        hmm.cutoffs.trusted = (cutoff, cutoff)
        loaded.append(hmm)
    
    Result = collections.namedtuple("Result", ["query", "hmm", "bitscore"])
    results = []
    for hits in pyhmmer.hmmsearch(loaded, proteins, bit_cutoffs="trusted"):
        hmm = hits.query_name.decode()
        for hit in hits:
            print(hit)
            if hit.included:
                results.append(Result(hit.name.decode(), hmm, hit.score))
    print(results)
    
    return results

def filter_hits(results, required):
    best_results = {}
    keep_query = set()
    for result in results:
        if result.query in best_results:
            previous_bitscore = best_results[result.query].bitscore
            if result.bitscore > previous_bitscore:
                best_results[result.query] = result
                keep_query.add(result.query)
            elif result.bitscore == previous_bitscore:
                if best_results[result.query].cog != hit.cog:
                    keep_query.remove(result.query)
        else:
            best_results[result.query] = result
            keep_query.add(result.query)
    hmm_hits = []
    for i in best_results:
        hmm_hits.append(best_results[i][1])
    if all(elem in hmm_hits for elem in required):
        for i in best_results:
            print(best_results[i][0], best_results[i][1], "{:.1f}".format(best_results[i][2]), sep="\t")
        return best_results
    else:
        print("Not all elements of the list are in the set.")

def parse_gbk(infile, parsed_df):
    """
    Function that retrieves information from a gbk file
    input:  annotated genome genbank file in GENBANK format;
            query, a string with protein accession
    output: returns organism, contig, and sequence information
    """

    #query = parsed_df['Accession'].iat[0]
    queries = parsed_df['Accession'].tolist()
    recs = [rec for rec in SeqIO.parse(infile, "genbank")]

    for i in parsed_df.index:
        query = parsed_df.at[i, 'Accession']
        for i in range(len(queries)):
            for rec in recs:
                for feat in rec.features:
                    if feat.type == "CDS":
                        try:       
                            if feat.qualifiers['locus_tag'][0] == queries[i]:
                                #print("Query is locus_tag")
                                nquery = query
                            elif feat.qualifiers['protein_id'][0] == queries[i]:
                                nquery = feat.qualifiers['locus_tag'][0]
                                #print("Query is protein_id")
                                #print("New query =",nquery)
                        except KeyError:
                            continue

            for rec in recs:
                feats_gene = [feat for feat in rec.features if feat.type == "gene"]
                feats_cds = [feat for feat in rec.features if feat.type == "CDS"]
                try:    
                    for feat in feats_cds:
                        locus_tag = feat.qualifiers['locus_tag'][0]       # get source (gene accession)
                        protein_id = feat.qualifiers['protein_id'][0]
                        if locus_tag == nquery:
                            organism = rec.annotations['organism']        # gets organism
                            contig = rec.id                               # contig is the rec.id
                            seq = feat.qualifiers['translation']
                            loc_cds = feat.location
                            #print(seq)
                            #print("Hit found:",locus_tag,"on contig",contig)
                            parsed_df.at[i, 'Organism'] = organism
                            parsed_df.at[i, 'Contig'] = contig
                            parsed_df.at[i, 'Sequence'] = str(seq)
                            parsed_df.at[i, 'location_cds'] = str(loc_cds)
                            parsed_df.at[i, 'Locus_tag'] = locus_tag
                            parsed_df.at[i, 'Protein_id'] = protein_id             
                except KeyError:
                    continue
                try:
                    for feat in feats_gene:            
                            locus_tag = feat.qualifiers['locus_tag'][0]
                            if locus_tag == nquery:
                                loc_mrna = feat.location
                                parsed_df.at[i, 'location_gene'] = str(loc_mrna)
                except KeyError:
                    continue
        
    return parsed_df

def update_cluster(cluster, best_results):
    """
    Adds predicted protein name to cluster Genbank files as 'product'
    """

    with open(cluster, "r") as handle:
        records = []
        for record in SeqIO.parse(handle, "genbank"):
            for feat in record.features:
                if feat.type == "CDS":
                    if feat.qualifiers["locus_tag"][0] in locus_tags:
                        feat.qualifiers['product'] = best_results[feat.qualifiers["locus_tag"][0]][1]
                        # Append the record to the list
                        records.append(record)
    SeqIO.write(records, cluster, "genbank") 

def get_range(infile, lt):
    """
    Takes locus_tags for two genes, gets range (the distance between these two genes), 
    and extracts 10 kb on either side of this range, outputting as a Genbank file
    Input: infile = target Genbank file;
           recs = parsed GBK file - output of parse_gbk;
           lt = locus_tag;
           n = iteratable integer for labelling purposes
    This code is adapted from "https://www.biostars.org/p/340270/" by user Joe, with some modifications
    """ 

    if not os.path.exists("clusters"):
        os.mkdir("clusters")

    file = str(infile.rsplit("/")[1].rsplit(".")[0])
    #print(file)
    # loop through GBK records and get information
    recs = [rec for rec in SeqIO.parse(infile, "genbank")]
    extracted_loci = []
    cluster_file = os.path.join("clusters", f"cluster_{file}_{str(len(lt))}.gb")
    for rec in recs:
        loci = [feat for feat in rec.features if feat.type == "CDS"]
        # try loop to catch exceptions if locus_tag doesn't exist
        try:
            start = min([int(l.location.start) for l in loci if 'locus_tag' in l.qualifiers and l.qualifiers['locus_tag'][0] in lt])
            end = max([int(l.location.start) for l in loci if 'locus_tag' in l.qualifiers and l.qualifiers['locus_tag'][0] in lt])
            (start and end)
            left_edge = start - 10000
            # can't have negative ints for start loc
            if left_edge < 0:
                left_edge = 0
            right_edge = end + 10000
            subrecord = []
            subrecord = rec[left_edge:right_edge]
            extracted_loci.append(subrecord)    
        except ValueError:
            continue    
        except NameError:
            print('Didn\'t get any indices even though the genes seemed to match. Couldn\'t slice.\n')
            pass
    #print(f"\tnumber of loci: {len(extracted_loci)}")
    if os.path.isfile(os.path.join("clusters", cluster_file)) is False:
        SeqIO.write(extracted_loci, cluster_file, "genbank")

def align_sequences(infile, out_filename):
    """
    Align sequences using MUSCLE
    """
    if not os.path.isfile(out_filename):
        cmd = f"muscle -in {infile} -out {out_filename}"
        subprocess.run([cmd], shell=True)
    return 

def calc_trusted_cutoffs(list_phmms, phmm_dir):
    """
    Calculates trusted cutoff thresholds for each pHMM.
    Input is a directory containing the pHMMs and a multi-Fasta file of the protein sequences used to generate the pHMMs.
    Returns dictionary with trusted cutoffs for each pHMM.
    -------
    """
    
    trusted_cutoffs = {}
    output_file = "trusted_cutoffs.txt"
    if os.path.isfile(output_file) is False or os.stat(output_file).st_size == 0:
        for hmm in list_phmms:
            prefix = hmm.split(".")[0]
            fasta_file = f"{prefix}.prots.fasta"
            results = multi_hmmsearch([hmm], fasta_file)                
            smallest_bitscore = min(results, key=lambda x: x.bitscore)
            print(smallest_bitscore[2])
                
            # calculate trusted cutoffs
            if len(results) >= 5:
                c = 0.2
            else:
                c = 0.5
            TC = smallest_bitscore[2] - (c * smallest_bitscore[2])
            print(TC)
            # write trusted cutoffs to dict and to file
            with open("trusted_cutoffs.txt", "a") as tc:
                tc.write(f"{prefix}\t{TC}\n")

    # will reuse existing trusted cutoff results file if it exists
    else:
        trusted_cutoffs = {}
        with open("trusted_cutoffs.txt", "r") as f:
            f_in = f.read()
            fs = f_in.split()
            for hmm in list_phmms:
                trusted_cutoffs["hmm"] = fs[0]
                trusted_cutoffs["TC"] = fs[1]
        print(trusted_cutoffs)

    if len(trusted_cutoffs) != len(list_phmms):
        print("Not all trusted cutoffs could be calculated for your pHMMs. Please check your inputs\n")
        exit(1)

    return trusted_cutoffs

##def search_gbk(i, gbk_file, hmm_list, phmm_dir, TC, required):
def search_gbk(i, gbk_file, phmms, required):
    """
    Overarching function that executes all the functions required to extract clusters
    from a GBK file using the input pHMM files.
    """
    
    gbk_prefix = gbk_file.split(".")[0].split("/")[1]
    cds_outfile = os.path.join("cds", f"{gbk_prefix}.fasta")
    gbk_results_file = f"{gbk_prefix}_results.tsv"
    global files_len
    get_cds(gbk_file, cds_outfile)
    print(f"Searching in {gbk_file:<30} - {i+1:5}/{files_len}")
    # if condition looks for existing TSV data file and reads data to df
    try:
        if os.path.exists(gbk_results_file) and os.stat(gbk_results_file).st_size > 0:
            #print("\tdata already written to tsv")
            df = pd.read_csv(gbk_results_file, sep='\t')
        else:
            # loop through each pHMM and search against GBK CDS
            hits = trusted_hmmsearch(phmms, cds_outfile, "trusted_cutoffs.txt")
            filtered_hits = filter_hits(hits, input_args.required)
            # then builds up contextual information for the hits from the GBK

            if len(filtered_hits) > 0:
                with open(gbk_results_file, "w") as f:
                    f.write(filtered_hits)

    except pd.errors.EmptyDataError:
        print(f"{gbk_results_file} is empty")
        return 

    if len(filtered_hits) > 0:
        locus_tags = df["query"].values.tolist()
        get_range(gbk_file, locus_tags)
    else:
        #print(f"\tNo clusters in {gbk_file}")
        return

    # Adds product names to each hit in each cluster genbank file
    clusters = os.listdir("clusters")
    n = 0
    for cluster in clusters:
        if cluster.split("_")[1] == gbk_prefix:
            #print(f"\tupdating cluster {cluster}")
            update_cluster(os.path.join("clusters", cluster), df)
            n += 1

def multiprocess_function(function, directory, var1, var2, var3, var4):
    """
    A function that uses the threading module to apply a 
    function to each file in a given directory
    """
    # Get a list of all files in the directory
    files = [f for f in os.listdir(directory) if f.endswith('.gbk') or f.endswith('.gb')]
    global files_len
    files_len = len(files)
    #print(files)
    # Determine the number of tasks to assign
    # Assigns four cpus per task
    if len(files)/4 < 1:
        cpus = int(multiprocessing.cpu_count()/len(files))
    else:
        cpus = int(multiprocessing.cpu_count()/4)
    #args_list = [(os.path.join(directory, file), var1, var2, var3, var4) for file in files]
    args_list = [(i, os.path.join(directory, file), var1, var2, var3, var4) for i, file in enumerate(files)]
    ##print(args_list)
    #print(f"{multiprocessing.cpu_count()} CPUs detected")
    #print(f"Split input GBKs into {batch_size} batches for parallelization")

    with multiprocessing.Pool(processes=cpus) as pool:
        # Apply a function to each file
        # using the worker processes in the pool, and
        # passing the additional input parameters to the
        # function
        #print("running pool")
        pool.starmap(function, args_list)

def main():

    print(f"{'':10}{'|-|-|-|- Cluster Search -|-|-|-|':^{30}}\n")
    args = get_args()
    phmm_dir = args.phmm_dir
    gbk_dir = args.gbk_dir
    try:
        required = (args.required).split(",")
        required_stripped = [os.path.splitext(x)[0] for x in required]
    except Exception:
        required_stripped = None

    if args.alignment is True:
        msas = [f for f in os.listdir(phmm_dir) if f.endswith('.aln')]
        for msa in msas:
            phmm_out = os.path.join(phmm_dir, f"{msa.split(".")[0]}.h3m")
            build_hmm(os.path.join(phmm_dir, msa), phmm_out)

    phmms = [f for f in os.listdir(phmm_dir) if f.endswith('.hmm') or f.endswith('.h3m')]
    #print(phmms)
    ##gbks = [f for f in os.listdir(gbk_dir) if f.endswith(('.gb', '.gbk'))]

    # calculate the trusted cutoff values for each pHMM
    ##trusted_cutoffs = calc_trusted_cutoffs(phmms, phmm_dir)

    if os.path.exists("hmm_output") is False:
        os.makedirs("hmm_output", exist_ok=True)
    ##for gbk in gbks:
    multiprocess_function(search_gbk, gbk_dir, phmms, required_stripped) ## phmm_dir, trusted_cutoffs, 

    clusters = os.listdir("clusters")
    # Reports the number of clusters extracted
    if len(clusters) == 1:
        print(f"\n{len(clusters)} cluster was found and is available in the \'clusters\' directory\n")
    elif len(clusters) == 0:
        print("\nNo clusters found :'(")
    else:
        print(f"\n{len(clusters)} clusters were found and are available in the \'clusters\' directory\n")

if __name__ == '__main__':
    main()

# Written by Kelly Styles
# Version 1