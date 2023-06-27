import os
import argparse
import subprocess
import sys
import pyhmmer
import multiprocessing
import collections
from Bio import SeqIO

version = 3
f"""
                    |-|-|-|- Cluster Search v{version} -|-|-|-|

This script will pull out sequence data containing hits to input pHMM models 
from a GBK file.  

Usage: 
    python3 cluster_search_v{str(version)}.py -g gbk_dir -p phmm_dir

Usage if you want to use a modifier for calculating very stringent trusted cutoffs: 
    python3 cluster_search_v{str(version)}.py -g gbk_dir -p phmm_dir -m 0.1

Usage if you want to use a modifier for calculating very relaxed trusted cutoffs: 
    python3 cluster_search_v{str(version)}.py -g gbk_dir -p phmm_dir -m 0.9

Usage if you want to filter hits to list of required pHMM hits: 
    python3 cluster_search_v{str(version)}.py -g gbk_dir -p phmm_dir -r list_required_models

Usage if you want to build pHMMs from multiple sequence alignments:
    python3 cluster_search_v{str(version)}.py -g gbk_dir -p phmm_dir -b

Usage if you want to align a multi-FASTA of proteins and build pHMMs from the multiple sequence alignments:
    python3 cluster_search_v{str(version)}.py -g gbk_dir -p phmm_dir -a -b

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
    parser.add_argument('-g', '--gbk_dir', action='store', type=str,
                        help='path to genbank file(s)', required=True)
    parser.add_argument('-p', '--phmm_dir', action='store', type=str,
                        help='path to pHMMs', required=True)
    parser.add_argument('-o', '--output_dir', action='store', type=str,
                        help='path to output directory', required=False)
    parser.add_argument('-b', '--build', action='store_true',
                        help='Build pHMM from multiple sequence alignments. Multiple sequence alignments must be in \'phmm_dir\' and have suffix \'.msa.fasta\'', required=False)
    parser.add_argument('-a', '--align', action='store_true',
                        help='Align multi-FASTA to make multiple sequence alignments. multi-FASTA must be in \'phmm_dir\' and have suffix \'.prots.fasta\'', required=False)
    parser.add_argument('-r', '--required', action='store', type=str,
                        help='comma-separated list of pHMMs required to be in the cluster, e.g., phmmA,phmmB,phmmC', required=False)
    parser.add_argument('-m', '--trusted_modifier', action='store', default=None, type=float,
                        help='modifier larger than 0 and including 1 for calculating trusted cutoff bitscores. 1 = no cutoff, 0.1 = very stringent. \
                        Default is 0.2 for alignments with >= 5 seqs and 0.5 for <= 4 seqs', required=False)
    parser.add_argument('-t', '--trusted_cutoffs_file', action='store', default=None, type=str,
                        help='file containing tab-separated trusted cutoffs for each pHMM, e.g. \'phmmA\tbitscore\'', required=False)
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        exit(1)
    
    return parser.parse_args()

def check_versions():
    """
    Print installed software versions
    """
    print("Installed Software:\n")

    # Python modules
    print(f"Python: {sys.version}")
    print(f"BioPython: {Bio.__version__}")
    
    # Bash modules
    cmd = "muscle --version | grep muscle | cut -d\"[\" -f1"
    subprocess.run(cmd, shell=True)

def get_cds(gbk, out_filename):
    """
    Input:  target GBK file
    Output: protein multi-FASTA file

    Extract protein sequences from a GBK file and write to a multi-FASTA file
    """
    if not os.path.isfile(out_filename):
        with open(gbk) as handle, open(out_filename, "w") as output_handle:
            for seq_record in SeqIO.parse(handle, "genbank"):
                for seq_feature in seq_record.features:
                    if seq_feature.type == "CDS":
                        try:
                            translation = seq_feature.qualifiers['translation'][0]
                            if len(translation) > 0:
                                # will take locus_tag preferentially, else
                                # will use protein_id
                                try:
                                    locus_tag = seq_feature.qualifiers['locus_tag'][0]
                                except KeyError:
                                    locus_tag = seq_feature.qualifiers['protein_id'][0]
                            output_handle.write(">%s \n%s\n" % (
                                locus_tag,     
                                translation))
                        except KeyError:
                            continue  # Skip the feature if 'translation' key is not present
    return out_filename

def build_hmm(msa, output_file):
    """
    Inputs: msa = multiple sequence alignment FASTA
    Output: output_file = pHMM file

    Builds a binary HMM 'h3m' file using pyhmmer from an alignment file
    """
    # load pyhmmer tools
    alphabet = pyhmmer.easel.Alphabet.amino()
    builder = pyhmmer.plan7.Builder(alphabet)
    background = pyhmmer.plan7.Background(alphabet)
    
    # load MSA
    with pyhmmer.easel.MSAFile(msa, digital=True, alphabet=alphabet) as msa_file:
        msa_d = msa_file.read()
    # add name to MSA
    prefix = msa.split("/")[-1].split(".")[0]
    msa_d.name = prefix.encode()
    # build pHMM
    hmm, _, _ = builder.build_msa(msa_d, background)

    # write pHMM to output_file
    with open(output_file, "wb") as f:
        print(f"Writing pHMM to {output_file}")
        hmm.write(f)
  
def multi_hmmsearch(hmm_list, fasta):
    """
    Inputs: hmms    = list of pHMMs
            fasta   = target protein multi-FASTA
    Output: results of hmmsearch

    Searches a list of pHMMs against a target protein multi-FASTA.
    """

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
    
    return results

def trusted_hmmsearch(hmms, fasta, trusted_cutoff_file):
    """
    Inputs: hmms    = list of pHMMs
            fasta   = target protein multi-FASTA
            trusted_cutoff_file = file with trusted cutoffs
    Output: results of hmmsearch

    Searches a list of pHMMs against a target protein multi-FASTA using trusted cutoffs.
    """
    alphabet = pyhmmer.easel.Alphabet.amino()
    # load trusted cutoffs
    cutoffs = {}
    with open(trusted_cutoff_file, "r") as tc_file:
        for line in tc_file:
            line = line.split()
            cutoffs[line[0]] = float(line[1])
    
    # load protein sequences
    with pyhmmer.easel.SequenceFile(fasta, digital=True, alphabet=alphabet) as seqs_file:
        proteins = seqs_file.read_block()
    #print(f"Loaded {len(proteins)} protein CDS from {fasta}")

    # load pHMMs and assign trusted cutoffs
    loaded = []
    for f in hmms:
        with pyhmmer.plan7.HMMFile(f) as hmm_file:
            hmm = hmm_file.read()
        try:
            cutoff = cutoffs[hmm.name.decode()]    
            hmm.cutoffs.trusted = (cutoff, cutoff)
            loaded.append(hmm)
        except KeyError:
            print("pHMMs names do not match the prefix, edit the NAME parameter in the pHMM file")
    
    # create collection to store results
    Result = collections.namedtuple("Result", ["query", "hmm", "bitscore"])
    
    # search loaded pHMMs against loaded proteins
    results = []
    for hits in pyhmmer.hmmsearch(loaded, proteins, bit_cutoffs="trusted"):
        hmm = hits.query_name.decode()
        for hit in hits:
            if hit.included:
                results.append(Result(hit.name.decode(), hmm, hit.score))
    
    return results

def filter_hits(results, gbk_results_file, required):
    """
    Input:  results = list of results
            gbk_results_file = path to output results file
            required = list of pHMMs that must be present in results
    Output: dictionary of filtered hits, also written to 'gbk_results_file'

    Will get the top hit for each pHMM, then check if all required pHMMS had hits.
    Writes these results to 'gbk_results_file'.
    """
    best_results = {}
    keep_query = set()

    # Filter out top hits for each pHMMs
    for result in results:
        if result.query in best_results:
            previous_bitscore = best_results[result.query].bitscore
            if result.bitscore > previous_bitscore:
                best_results[result.query] = result
                keep_query.add(result.query)
            elif result.bitscore == previous_bitscore:
                keep_query.remove(result.query)
        else:
            best_results[result.query] = result
            keep_query.add(result.query)
    hmm_hits = [best_results[i][1] for i in best_results]

    # Check all require pHMMs have hits
    if len(hmm_hits) > 0:
        if required is not None:
            if all(elem in hmm_hits for elem in required):
                with open(gbk_results_file, "w") as f:
                    for i in best_results:
                        line = "\t".join([str(best_results[i][0]), str(best_results[i][1]), "{:.1f}".format(best_results[i][2])])
                        f.write(line + "\n")
                return best_results
        else:
            with open(gbk_results_file, "w") as f:
                for i in best_results:
                    line = "\t".join([str(best_results[i][0]), str(best_results[i][1]), "{:.1f}".format(best_results[i][2])])
                    f.write(line + "\n")
            return best_results

def open_gbk(infile):
    """
    Input:  target GBK file
    Output: parsed GBK object

    Parses GBK usng SeqIO from BioPython.
    """
    recs = [rec for rec in SeqIO.parse(infile, "genbank")]
    return recs

def get_range(infile, recs, lt):
    """
    Input: infile = target Genbank file
           recs = parsed GBK file - output of open_gbk
           lt = list of locus tags
    Output: output GBK file containing contextual hits of locus tags

    Takes locus_tags for two genes, gets range (the distance between these two genes), 
    and extracts 10 kb on either side of this range, outputting as a Genbank file
    This function has been adapted from "https://www.biostars.org/p/340270/" by user Joe.
    """ 

    file = str(infile.rsplit("/")[1].rsplit(".")[0])
    #print(file)
    # loop through GBK records and get information

    extracted_loci = []
    cluster_file = os.path.join(out_dir, "clusters", f"cluster_{file}_{str(len(lt))}.gb")
    for rec in recs:
        loci = [feat for feat in rec.features if feat.type == "CDS"]
        organism = rec.annotations['organism']
        # try loop to catch exceptions if locus_tag doesn't exist
        try:
            start = min([int(l.location.start) for l in loci if l.qualifiers['locus_tag'][0] in lt])
            end = max([int(l.location.end) for l in loci if l.qualifiers['locus_tag'][0] in lt])
            (start and end)
            #print(rec.id, start, end)
            if start - 10000 < 0:
                left_edge = 0
            else:
                left_edge = start - 10000
            right_edge = end + 10000
            #print(rec.id, left_edge, right_edge)
            subrecord = []
            subrecord = rec[left_edge:right_edge]
            organism_ = organism.replace(" ","_")
            extracted_loci.append(subrecord)   
            #print(len(extracted_loci)) 
        except KeyError:
            continue
        except ValueError:
            continue    
        except NameError:
            #print('Didn\'t get any indices even though the genes seemed to match. Couldn\'t slice.\n')
            pass
    #print(f"number of loci: {len(extracted_loci)}")
    if os.path.isfile(os.path.join(out_dir, "clusters", cluster_file)) is False:
        SeqIO.write(extracted_loci, cluster_file, "genbank")

def align_sequences(infile, out_filename):
    """
    Input:  protein multi-FASTA
    Output: alignment FASTA file

    Align sequences using MUSCLE
    """
    if not os.path.isfile(out_filename):
        print(f"aligning {out_filename}")
        cmd = f"muscle -in {infile} -out {out_filename}"
        subprocess.run([cmd], shell=True)
    return 

def calc_trusted_cutoffs(list_phmms, phmm_dir, trusted_modifier, trusted_cutoff_input):
    """
    Calculates trusted cutoff thresholds for each pHMM.
    Input is a directory containing the pHMMs and a multi-Fasta file of the protein sequences used to generate the pHMMs.
    Returns dictionary with trusted cutoffs for each pHMM.
    -------
    """
    
    trusted_cutoffs = []
    # if user supplied trusted cutoffs exists, use instead
    if trusted_cutoff_input is not None:
        output_file = trusted_cutoff_input
    else:
        output_file = "trusted_cutoffs.txt"
    
    # Calculate trusted cutoffs for each pHMM,
    # else reuse existing trusted cutoff results file 'output_file' if it exists
    if os.path.isfile(output_file) is False or os.stat(output_file).st_size == 0 or trusted_modifier is not None:
        print("Calculating the trusted cutoffs for each HMM...")
        for hmm in list_phmms:
            tc_dict = {}
            prefix = hmm.split("/")[-1].split(".")[0]
            fasta_file = os.path.join(phmm_dir, f"{prefix}.prots.fasta")
            hmm_path = os.path.join(phmm_dir, hmm)
            # search pHMM against the multi FASTA used to generate the pHMM
            results = multi_hmmsearch([hmm_path], fasta_file)              
            # get lowest scoring hit  
            smallest_bitscore = min(results, key=lambda x: x.bitscore)
                
            # Calculate trusted cutoffs using modifier
            # Modifier can come from input args, or use defaults from 
            # Villebro et al., 2019 paper, based on number of seqs in alignment
            # The larger the 'c' value, the less stringent the trusted cutoff
            if trusted_modifier is not None:
                m = trusted_modifier
            elif len(results) >= 5:
                m = 0.2
            else:
                m = 0.5
            trusted = smallest_bitscore[2] - (m * smallest_bitscore[2])
            
            # write trusted cutoffs to dict and to file
            tc_dict["pHMM"] = smallest_bitscore[1]
            tc_dict["trusted"] = trusted  
            with open("trusted_cutoffs.txt", "a") as tc:
                tc.write(f"{prefix}\t{trusted}\n")
    else:
        tc_dict = {}
        print("* Reusing previously calculated trusted cutoffs for each pHMM.")
        with open("trusted_cutoffs.txt", "r") as f:
            lines = f.readlines()
            for line in lines:
                line = line.split()
                for phmm in list_phmms:
                    if line[0] == phmm.split(".")[0]:
                        tc_dict["pHMM"] = line[0]
                        tc_dict["trusted"] = line[1]
                        trusted_cutoffs.append(tc_dict)

    if len(trusted_cutoffs) != len(list_phmms):
        print(f"TCs:{len(trusted_cutoffs)} vs List pHMMs:{len(list_phmms)}")
        print("Not all trusted cutoffs could be calculated for your pHMMs. Please check your inputs\n")

    return trusted_cutoffs

def update_cluster(cluster, best_hits_dict):
    """
    Adds predicted protein name to cluster Genabnk files as 'product'
    """
    # get locus tags to lookup
    locus_tags = [best_hits_dict[i][0] for i in best_hits_dict]

    # open the cluster GBK file, loop through
    with open(cluster, "r") as handle:
        records = []
        for record in SeqIO.parse(handle, "genbank"):
            for feat in record.features:
                if feat.type == "CDS":
                    # if locus tag matches one in 'locus_tags' add phmm used and bitscore to the record
                    if feat.qualifiers["locus_tag"][0] in locus_tags:
                        feat.qualifiers['product'] = best_hits_dict[feat.qualifiers["locus_tag"][0]][1]
                        feat.qualifiers['bitscore'] = best_hits_dict[feat.qualifiers["locus_tag"][0]][2]
            # Append the updated record to the list
            records.append(record)
    SeqIO.write(records, cluster, "genbank") 
   
def search_gbk(i, gbk_file, phmms, required, phmm_dir):
    """
    Overarching function that executes all the functions required to extract clusters
    from a GBK file using the input pHMM files.
    """
    
    gbk_prefix = gbk_file.split(".")[0].split("/")[-1]
    cds_outfile = os.path.join(out_dir, "cds", f"{gbk_prefix}.fasta")
    gbk_results_file = os.path.join(out_dir, "hmm_output", f"{gbk_prefix}_results.tsv")
    global files_len
    hits = None
    get_cds(gbk_file, cds_outfile)
    gbk_file_formatted = gbk_file[:35].ljust(35+5)
    print(f" - {gbk_file_formatted} - {i+1:5}/{files_len}")
    
    # if condition looks for existing TSV data file and reads data to df
    try:
        if os.path.exists(gbk_results_file) and os.stat(gbk_results_file).st_size > 0:
            hits = []
            with open(gbk_results_file) as f:
                lines = f.readlines()
                for line in lines:
                    #print(line)
                    line = line.split()
                    hits.append(line[0], line[1], line[2])
        else:
            phmm_paths = []
            for phmm in phmms:
                path = os.path.join(phmm_dir, phmm)
                phmm_paths.append(path)
            # loop through each pHMM and search against GBK CDS
            hits = trusted_hmmsearch(phmm_paths, cds_outfile, "trusted_cutoffs.txt")
    except TypeError:
        return 

    # extracts hit loci with 10 kb context from GBK file
    #print(hits)
    if hits is not None:
        # open GBK records
        recs = open_gbk(gbk_file)
        # filter hits
        filtered_hits = filter_hits(hits, gbk_results_file, required)
        #print(filtered_hits)
        if filtered_hits is not None:
            locus_tags = [filtered_hits[i][0] for i in filtered_hits]
            if len(locus_tags) > 0:
                #print(f"extracting cluster for {gbk_file}")
                get_range(gbk_file, recs, locus_tags)

        # Adds product names to each hit in each cluster genbank file
        clusters = os.listdir(os.path.join(out_dir, "clusters"))
        for cluster in clusters:
            if cluster.split("_")[1] == gbk_prefix:
                update_cluster(os.path.join(out_dir, "clusters", cluster), filtered_hits)

def multiprocess_function(function, directory, var1, var2, var3):#, var4):
    """
    A function that uses the threading module to apply a   
    function to each file in a given directory
    """
    # Get a list of all files in the directory
    files = [f for f in os.listdir(directory) if f.split(".")[-1].startswith('gb')]
    global files_len
    files_len = len(files)

    # Determine the number of tasks to assign
    # Assigns four cpus per task
    if len(files)/4 < 1:
        cpus = int(multiprocessing.cpu_count()/len(files))
    else:
        cpus = int(multiprocessing.cpu_count()/4)
    print(f"{GREEN}* {multiprocessing.cpu_count()} CPUs detected {END}")
    print(f"{GREEN}* Split input GBKs into {cpus} batches for parallelization {END}\n")
    print(f"Searching for hits in {len(files)} files:")

    # Map arugments to apply to function
    args_list = [(i, os.path.join(directory, file), var1, var2, var3) for i, file in enumerate(files)]
    with multiprocessing.Pool(processes=cpus) as pool:
        # Apply a function to each file
        # using the worker processes in the pool, and
        # passing the additional input parameters to the
        # function
        #print("running pool")
        pool.starmap(function, args_list)

def make_path(path):
    """
    Makes a path, regardless of if the parent dir exists or not.
    Functionally equivalent to `mkdir -p` in BASH
    """
    try:
        os.makedirs(path)
    except FileExistsError:
        pass

# Some functions that print pretty terminal output using ANSI codes
GREEN = "\u001b[32m"
BLUE = "\u001b[34m"
RED = "\033[91m"
YELLOW = "\u001b[33m"
BOLD = "\033[1m"
END = "\033[00m"

def main():

    print(f"\n{BLUE} {BOLD} {'':10}{'|-|-|-|- Cluster Search -|-|-|-|':^{30}} {END}\n")
    #check_versions()
    args = get_args()
    phmm_dir = args.phmm_dir
    gbk_dir = args.gbk_dir
    global out_dir
    out_dir = args.output_dir
    trusted_mod = args.trusted_modifier
    trusted_file = args.trusted_cutoffs_file
    try:
        required = (args.required).split(",")
        required_stripped = [os.path.splitext(x)[0] for x in required]
    except Exception:
        required_stripped = None

    if args.align is True:
        multis = [f for f in os.listdir(phmm_dir) if f.endswith('.prots.fasta')]
        for multi in multis:
            prefix = multi.split("/")[-1].split(".")[0]
            multi_out = os.path.join(phmm_dir, f"{prefix}.msa.fasta")
            align_sequences(os.path.join(phmm_dir, multi), multi_out)

    if args.build is True:
        msas = [f for f in os.listdir(phmm_dir) if f.endswith('.msa.fasta')]
        for msa in msas:
            prefix = msa.split("/")[-1].split(".")[0]
            phmm_out = os.path.join(phmm_dir, f"{prefix}.h3m")
            build_hmm(os.path.join(phmm_dir, msa), phmm_out)
    
    make_path(os.path.join(out_dir, "cds"))
    make_path(os.path.join(out_dir, "hmm_output"))
    make_path(os.path.join(out_dir, "clusters"))

    phmms = [f for f in os.listdir(phmm_dir) if f.endswith('.hmm') or f.endswith('.h3m')]

    # calculate the trusted cutoff values for each pHMM
    calc_trusted_cutoffs(phmms, phmm_dir, trusted_mod, trusted_file)
    
    # apply sarch_gbk function to each GBK file
    multiprocess_function(search_gbk, gbk_dir, phmms, required_stripped, phmm_dir)#, trusted_cutoffs, 

    clusters = os.listdir(os.path.join(out_dir, "clusters"))
    # Reports the number of clusters extracted 
    if len(clusters) == 1:
        print(f"\n{len(clusters)} cluster was found and is available in the \'clusters\' directory\n")
    elif len(clusters) == 0:
        print(f"\n{RED} {BOLD} No clusters found :'( {END}")
    else:
        cluster_path = os.path.join(out_dir, "clusters")
        print(f"\n{YELLOW} {BOLD} {len(clusters)} clusters were found and are available in the \'{cluster_path}\' directory {END}\n")

if __name__ == '__main__':
    main()

# Written by Kelly Styles
# Version 3