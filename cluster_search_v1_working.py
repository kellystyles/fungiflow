import os
import argparse
import subprocess
import sys
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

def make_blastdb(infile):
    """
    Create a BLAST database from a multi-FASTA protein sequence file
    """
    db_name = infile.rsplit(".", 1)[0].split("/")[1]
    if os.path.exists("blast_db") == False:
        os.mkdir("blast_db")
    
    if not os.path.isfile(f"blast_db/{db_name}.pin"):
        cmd = f"makeblastdb -in {infile} -input_type fasta -dbtype prot -out blast_db/{db_name}"
        subprocess.run([cmd], shell=True)
    return f"blast_db/{db_name}"

def blastp(query, db, out_filename):
    """
    Run BLASTp on a protein sequence file against a BLAST database
    """
    if not os.path.isfile(out_filename):
        cmd = f"blastp -query {query} -db {db} -outfmt \"6 std\" -out {out_filename}"
        subprocess.run([cmd], shell=True)
    return out_filename

def hmmsearch(hmm, fasta, out_filename):
    """
    Run HMMsearch on a HMM file against a multi-FASTA protein sequence database
    """
    if os.path.isfile(out_filename) is False:
        cmd = ["hmmsearch", "-E", "1e-5", "--noali", "--notextw", "--cpu","4", "-o", out_filename, hmm, fasta]
        subprocess.run(cmd)

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

def update_cluster(cluster, sorted_df):
    """
    Adds predicted protein name to cluster Genabnk files as 'product'
    """
    locus_tags = sorted_df["Locus_tag"].to_list()
    with open(cluster, "r") as handle:
        records = []
        for record in SeqIO.parse(handle, "genbank"):
            for feat in record.features:
                if feat.type == "CDS":
                    if feat.qualifiers["locus_tag"][0] in locus_tags:
                        feat.qualifiers['product'] = sorted_df.at[sorted_df[sorted_df['Locus_tag'] == feat.qualifiers['locus_tag'][0]].index[0], 'phmm']
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
#    files = os.listdir(list_phmms)
    output_file = "trusted_cutoffs.txt"
    #print(list_phmms)
    if os.path.isfile(output_file) is False or os.stat(output_file).st_size == 0:
        for file in list_phmms:
            #print(file)
            results = pd.DataFrame()
            clean_coll = []
            prefix = file.split(".")[0]
            fasta_file = os.path.join("hmms", f"{prefix}.fasta")
            results_file = os.path.join("hmms", f"{prefix}.txt")

            # hmmsearch the input pHMM file against the input FASTA file
            hmmsearch(os.path.join(phmm_dir, file), fasta_file, results_file)
            # read the output file from above and obtain hit lines
            forbidden = ['!', '?', '*', 'PP']
            with open(results_file, 'r') as read_obj:
                for line in read_obj:
                    cleaned = line.strip()
                    if cleaned.startswith(tuple('0123456789')) and all(forbidden_char not in cleaned for forbidden_char in forbidden):
                        clean_coll.append(cleaned.split(maxsplit=9)[:9])
                
            # add data to pandas dataframe
            results = pd.DataFrame(clean_coll)
            results.columns = ['seq_E-value', 'seq_score', 'seq_bias', 'dom_E-value', 'dom_score', 'dom_bias', 'exp', 'N', 'Accession']   
                
            seq_min = float(results.seq_score.min())
            dom_min = float(results.dom_score.min())
                
            # calculate trusted cutoffs
            if len(results) >= 5:
                c = 0.2
            else:
                c = 0.5
            seq_TC = seq_min - (c * seq_min)
            dom_TC = dom_min - (c * dom_min)
                
            # write trusted cutoffs to dict and to file
            trusted_cutoffs[prefix]={'Name':prefix,'seq_TC':round(seq_TC, 2),'dom_TC':round(dom_TC, 2)}
            with open("trusted_cutoffs.txt", "w") as tc:
                tc.write(str(trusted_cutoffs))

    # will reuse existing trusted cutoff results file if it exists
    else:
        trusted_cutoffs = {}
        with open("trusted_cutoffs.txt", "r") as f:
            f_in = f.read()
            trusted_cutoffs = ast.literal_eval(f_in)

    if len(trusted_cutoffs) != len(list_phmms):
        #print("Not all trusted cutoffs could be calculated for your pHMMs. Please check your inputs\n")
        exit(1)

    return trusted_cutoffs

def search_output(search_results, gbk_file, phmm, seq_TC, dom_TC):
    """
    Filters out hits below the trusted cutoff bitscores, then outputs results to a dictionary.
    
    Parameters:
        search_results (str): HMMSEARCH output file.
        seq_TC (float): Trusted cutoff bitscore for sequence.
        dom_TC (float): Trusted cutoff bitscore for domain.
    
    Returns:
        parsed (dict): A nested dictionary with hit information, bitscores, and e-values for sequence and domain hits.
    """
    df = pd.DataFrame()
    clean_coll = []
    multiplier = 0.5

    # Read in input file and obtain hit lines.
    with open(search_results, 'r') as read_obj:
        for line in read_obj:
            cleaned = line.strip()
            if (
                cleaned.startswith(tuple('0123456789'))
                and ('!' not in cleaned)
                and ('?' not in cleaned)
                and ('*' not in cleaned)
                and ('PP' not in cleaned)
            ):
                clean_coll.append(cleaned.split(maxsplit=9)[:9])
            elif cleaned == "[No hits detected that satisfy reporting thresholds]":
                hmm = str(search_results).split("_")[0]
                #print(f"No hits detected that satisfy reporting thresholds were detected for {hmm}")
    try:
        df = pd.DataFrame(clean_coll)
        df.columns = ['seq_E-value', 'seq_score', 'seq_bias', 'dom_E-value', 'dom_score', 'dom_bias', 'exp', 'N', 'Accession']
        df["input"] = gbk_file
        df["phmm"] = phmm.split(".")[0]

        # Filter out hits that have lower sequence and domain bitscores than the TCs.
        df = df[df['seq_score'].astype(float) >= seq_TC * multiplier]
        df = df[df['dom_score'].astype(float) >= dom_TC * multiplier]
        #print(df)
    except ValueError:
        #print(f"No hits for {phmm} in {gbk_file}")
        pass
    if len(df) > 0:
        return df

def search_gbk(i, gbk_file, hmm_list, phmm_dir, TC, required):
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
            # then builds up contextual information for the hits from the GBK
            # writes df to a TSV data file
            new_df = pd.DataFrame()
            for phmm in hmm_list:
                print(phmm)
                hmm_prefix = phmm.split(".")[0]
                search_outfile = os.path.join("hmm_output", f"{hmm_prefix}_{gbk_prefix}.txt")
                if os.path.isfile(search_outfile) is False:
                    hmmsearch(os.path.join(phmm_dir, phmm), cds_outfile, search_outfile)
                parsed = search_output(search_outfile, gbk_file, phmm, TC[hmm_prefix]['seq_TC'], TC[hmm_prefix]['dom_TC'])
                
                print(parsed)
                if parsed is not None:
                    print("LEn of parsed", len(parsed))
                    parse_gbk(gbk_file, parsed)
                    new_df = pd.concat([new_df, parsed])
            if len(new_df) > 0:
                new_df.to_csv(gbk_results_file, index=False, sep="\t")
                df = pd.read_csv(gbk_results_file, sep='\t')
            else:
                #print(f"\tNo pHMM hits in {gbk_file}")
                return
    except pd.errors.EmptyDataError:
        print(f"{gbk_results_file} is empty")
        return 
    # prepares a new df containing rows with max scores for each pHMM
    #print(df)
    if df["seq_score"].dtypes != 'float64':
        df["seq_score"] = pd.to_numeric(df["seq_score"], errors='coerce')
    sorted_df = df.loc[df.groupby("phmm")["seq_score"].idxmax()]

    # extracts each locus_tag and its context from GBK infile  
    # writes data to a cluster GBK file
    phmm_vals = sorted_df["phmm"].values.tolist()
    phmm_vals_stripped = [os.path.splitext(x)[0] for x in phmm_vals]
    locus_tags = sorted_df["Locus_tag"].values.tolist()
    if required is not None:
        #print(phmm_vals_stripped, required)
        result = all(x in phmm_vals_stripped for x in required)
    else:
        result = True
    print(result)
    if result == True:
        print(locus_tags)
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
            update_cluster(os.path.join("clusters", cluster), sorted_df)
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
    phmms = [f for f in os.listdir(phmm_dir) if f.endswith('.hmm')]
    #print(phmms)
    ##gbks = [f for f in os.listdir(gbk_dir) if f.endswith(('.gb', '.gbk'))]

    # calculate the trusted cutoff values for each pHMM
    trusted_cutoffs = calc_trusted_cutoffs(phmms, phmm_dir)

    if os.path.exists("hmm_output") is False:
        os.makedirs("hmm_output", exist_ok=True)
    ##for gbk in gbks:
    multiprocess_function(search_gbk, gbk_dir, phmms, phmm_dir, trusted_cutoffs, required_stripped)

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