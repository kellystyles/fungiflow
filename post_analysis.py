import os
import sys
import subprocess
import datetime
import shutil
import json
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import library as lib

sns.set_palette('tab10')
sns.set(rc={'figure.figsize':(12,10)})
sns.set_style(style='white')
sns.set_context("paper", font_scale=1)

def itsx(input_args,filenames,itsx_path):
    """
    Extracts the ITS sequences from an assembly FASTA file using `itsx`.
    
    Input:      assembly FASTA file
    Output:     itsx directory with ITS FASTA files
    """
    stdout = os.path.join(itsx_path,f"{input_args.array}.out")
    stderr = os.path.join(itsx_path,f"{input_args.array}.err")

    cmd = ["ITSx","-i",filenames.assembly_fasta,"-o",os.path.join(itsx_path,input_args.array),"-t","F","--cpu",input_args.cpus,"--preserve","--only_full","T","--temp",itsx_path]
    if len(filenames.singularity) > 0: cmd = filenames.singularity + cmd
    if len(input_args.nanopore) > 0: cmd = cmd + ["--nhmmer", "T"]
    try:
        lib.print_n("Extracting the ITS sequence with ITSx")
        lib.execute(cmd,stdout,stderr)
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output)  

def blastn(input_args,filenames,blast_path):
    """
    Will run `blastn` on an ITS FASTA file and output to a results TXT file.

    Input:      ITS input FASTA file
    Output:     blastn TXT output file    
    """

    stdout = os.path.join(blast_path,f"{input_args.array}.out")
    stderr = os.path.join(blast_path,f"{input_args.array}.err")

    cmd = f"blastn -db {input_args.its_db} -query {filenames.its_fasta} -num_threads {input_args.cpus} -outfmt \"6 qacc sscinames sacc pident bitscore evalue stitle\" -max_target_seqs 5 -out {filenames.blast_out}"
    if len(filenames.singularity) > 0: cmd = str(" ".join(filenames.singularity)) + " " + cmd
    commands = [cmd]
    try:
        lib.print_n(f"BLASTn of {filenames.its_fasta} against {input_args.its_db}")
        lib.execute_shell(commands,stdout,stderr)
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output)

def parse_blast(results_path,array):
    """
    Parses BLAST output files and appends the extracted information to a pandas DataFrame.

    Input:      BLAST output file from blast[x]() function [path]
    Output:     pandas DataFrame object
    """

    columns = ["Query","Species","Accession","Identity","Bitscore","E-value","Subject Name"]

    if os.path.exists(results_path):
        lib.print_n(f"Parsing BLAST results for {results_path}")
        try:
            dataframe = pd.read_csv(results_path, sep="\t", header=None)
            dataframe.columns = columns
            dataframe = dataframe.head(1)
            dataframe["ID"] = array
            print(dataframe)
            return dataframe
        except pd.errors.EmptyDataError:
            print(f"No BLAST results obtained for {results_path}")
    else:
        print(f"The BLAST for {results_path} was not completed successfully")
        
def quast(input_args,filenames,output_path):
    """
    Runs `quast` genome analysis software on an assembly.

    Input:      assembly FASTA file and GFF file (optional)
    Output:     QUAST output folder
    """

    stdout = os.path.join(output_path,f"{input_args.array}.out")
    stderr = os.path.join(output_path,f"{input_args.array}.err")

    try:
        if os.path.exists(filenames.funannotate_gff) == True and os.stat(filenames.funannotate_gff).st_size > 0:
            lib.print_n(f"Funannotate annotations GFF file exists, analysing {filenames.funannotate_gff} and {filenames.assembly_fasta} with quast")
            cmd = ["quast.py",filenames.assembly_fasta,"-o",output_path,"-g",filenames.funannotate_gff,"-t",input_args.cpus,"--fungus","-L","-b","--pe1",filenames.trimmedf,"--pe2",filenames.trimmedr]
            if len(filenames.singularity) > 0: cmd = filenames.singularity + cmd
        else:
            lib.print_n(f"Funannotate annotations GFF file does not exist, analysing {filenames.assembly_fasta} with quast")
            cmd = ["quast.py",filenames.assembly_fasta,"-o",output_path,"-t",input_args.cpus,"--fungus","-L","-b","--pe1",filenames.trimmedf,"--pe2",filenames.trimmedr]
            if len(filenames.singularity) > 0: cmd = filenames.singularity + cmd
    except AttributeError:
        lib.print_n(f"Funannotate annotations GFF file does not exist, analysing {filenames.assembly_fasta} with quast")
        cmd = ["quast.py",filenames.assembly_fasta,"-o",output_path,"-t",input_args.cpus,"--fungus","-L","-b","--pe1",filenames.trimmedf,"--pe2",filenames.trimmedr]
        if len(filenames.singularity) > 0: cmd = filenames.singularity + cmd

    try:
        lib.execute(cmd,stdout,stderr)
        lib.file_exists_exit(filenames.quast_report,"Quast successfully analysed the assembly!","Quast failed... check the logs and your inputs")
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output) 

def parse_quast(filenames,array):
    """
    Parses `quast` output reports, returning a pandas dataframe.

    Input:      quast TSV report
    Output:     pandas DataFrame object
    """

    lib.print_n(f"Parsing Quast results for {filenames.quast_report}")
    try:
        dataframe = pd.read_csv(filenames.quast_report, sep="\t", header=0)
        dataframe["ID"] = array
        print(dataframe)
        return dataframe
    except pd.errors.EmptyDataError:
        lib.print_e(f"The QUAST report {filenames.quast_report} was not parsed successfully")

def antismash(input_args,filenames,antismash_out):
    """
    Runs `antismash` BGC analysis on assembly file.
        - antismash will have a cry about files already being present in the output folder, 
          so STDOUT and STDERR files found in parent dir, rather than 'antismash_out'.
    
    Input:      annotated assembly GBK file
    Output:     antismash output folder
    """
    
    stdout = os.path.join(input_args.directory_new,f"{input_args.array}_antismash.out")
    stderr = os.path.join(input_args.directory_new,f"{input_args.array}_antismash.err")
    
    cmd = ["antismash","--cpus",input_args.cpus,"--taxon","fungi","--fullhmmer","--genefinding-tool","glimmerhmm","--cc-mibig","--smcog-trees","--cb-general","--cb-subclusters","--cb-knownclusters","--minlength","5000","--output-dir",antismash_out,"--output-basename",input_args.array,filenames.antismash_assembly]
    if len(filenames.antismash) > 0: cmd = filenames.antismash + cmd
    try:
        lib.print_n(f"Predicting gene clusters from {filenames.antismash_assembly} with antiSMASH")
        lib.execute(cmd,stdout,stderr)
        lib.file_exists(filenames.antismash_json,"antiSMASH successfully analysed the assembly!","antiSMASH failed... check the logs and your inputs")
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output)

def open_json(json_file):
    """
    Opens a JSON file to an object
    Input: JSON file
    Output: JSON object
    """
    try:
        with open(json_file, 'r') as f:
            try:
                contents = json.load(f)
                return contents
            except json.JSONDecodeError:
                print(f"JSON file {json_file} does not exist")
                return None
    except:
        print(f"error occurred while parsing {json_file}:", sys.exc_info()[1])
        return None

def get_location(location):
    """
    Extracts coordinates of a BGC from an opened antiSMASH JSON record
    Input: JSON object location record
    Output: location coordinates
    """
    loc_ = location.replace("<", "").replace(">", "").split(":")
    # adds 1 to the starting base value to be inclusive of that base
    return int(loc_[0].lstrip("[")), int(loc_[1].rstrip("]"))

def plot_violin(dataframe, input_args, results_path):
    """
    Plots some figures summarizing the BGC information from parsed antiSMASH
    data using the seaborn package and save it as an SVG file:
     - violinplot of each BGC category identified vs. BGC length
    """
    
    ax = sns.violinplot(data=dataframe, x=dataframe["BGC_length"], y=dataframe["BGC_type"], group_by=dataframe["BGC_type"], palette="tab10", dodge=False)
    ax.get_figure().savefig(os.path.join(results_path,f"{input_args.array}_bgctype_vs_bgclength.svg"))
    plt.clf()

def plot_strip(dataframe, input_args, results_path):
    """
    Plots some figures summarizing the BGC information from parsed antiSMASH
    data using the seaborn package and save it as an SVG file:
     - swarmplot of contig edge (True/False) vs BGC length, with known clusters 
       labelled
    """
    ax = sns.swarmplot(data=dataframe, x=dataframe["BGC_length"], y=dataframe["partial"], hue=dataframe["BGC_type"], palette="tab10")
    ax.legend(loc='right', bbox_to_anchor=(1.5, 0.5), title="BGC type")
    flag = True
    for line in range(0,dataframe.shape[0]):
        if flag == True:
            align = "left"
            offset = -0.02
            flag = False
        else:
            align = "right"
            offset = 0.01
            flag = True
        ax.text(x=dataframe.BGC_length[line], y=offset, s=str(dataframe.kc_desc[line].split("/")[0]), \
                horizontalalignment=align, color='black', rotation=45, \
                yourfontsize=8, rotation_mode="anchor")
    ax.get_figure().savefig(os.path.join(results_path,f"{input_args.array}_contigedge_vs_bgclength_kc.svg"))
    plt.clf()

def parse_json(json_object,json_file):
    """
    Extracts information about the BGCs from a JSON object,
    including knownclusterblast data.
    Input: JSON object (from open_json), JSON file (str)
    Output: pandas array of information
    """
    bgc_list = []
    array = pd.DataFrame()
    bgc_num = 1
    types_c = 0
    input_file = str(json_file.split("\\")[-1].split(".")[0])
    for i in range(len(json_object['records'])):
        contig = json_object['records'][i]['id'].split("_")[0]
        # obtains a list of the BGC types searched for in the records
        if "antismash.detection.hmm_detection" in json_object['records'][i]["modules"]:
            while types_c < 1:
                types = json_object['records'][i]["modules"]["antismash.detection.hmm_detection"]["enabled_types"]
                #print(types)
                types_c += 1
        # loops through candidate clusters and extract information about each cluster
        for j in range(len(json_object['records'][i]['features'])):
            record = json_object['records'][i]['features'][j]
            kc_record = json_object['records'][i]["modules"]
            if 'cand_cluster' in record['type']:
                # prepares dict and adds information to this dict
                a = {}
                product = record["qualifiers"]['product']
                edge = record["qualifiers"]["contig_edge"][0]
                bgc_begin, bgc_end = get_location(record["location"])
                bgc_length = bgc_end - bgc_begin
                a["file"] = input_file
                a["contig"] = contig
                a["partial"] = edge
                a["BGC_type"] = product
                a["BGC_location"] = record["location"].lstrip("[").rstrip("]")#f"{bgc_begin}:{bgc_end}"
                a["BGC_length"] = bgc_length
                a["protocluster"] = record["qualifiers"]["protoclusters"][0]
                a["kind"] = record["qualifiers"]["kind"][0]
                # extracts the known cluster information
                # adapted from the above notebook
                if "antismash.modules.clusterblast" in kc_record:
                    for kc in kc_record["antismash.modules.clusterblast"]["knowncluster"]["results"]:
                        kc_total_hits = int(kc["total_hits"])
                        if kc_total_hits > 0:
                            a["kc_total_hits"] = kc_total_hits
                            a["kc_mibig_acc"] = kc["ranking"][0][0]["accession"]
                            a["kc_desc"] = kc["ranking"][0][0]["description"]
                            a["kc_type"] = kc["ranking"][0][0]["cluster_type"]
                            a["kc_hits"] = int(kc["ranking"][0][1]["hits"])
                            a["kc_core_gene_hits"] = int(kc["ranking"][0][1]["core_gene_hits"])
                            a["kc_blast_score"] = int(kc["ranking"][0][1]["blast_score"])
                            a["kc_synteny_score"] = int(kc["ranking"][0][1]["synteny_score"])
                            a["kc_core_bonus"] = int(kc["ranking"][0][1]["core_bonus"])
                        else:
                            a["kc_total_hits"] = kc_total_hits
                            a["kc_mibig_acc"] = ""
                            a["kc_desc"] = ""
                            a["kc_type"] = ""
                            a["kc_hits"] = 0
                            a["kc_core_gene_hits"] = 0
                            a["kc_blast_score"] = 0
                            a["kc_synteny_score"] = 0
                            a["kc_core_bonus"] = 0                           
                    bgc_list.append(a)
                    # might need to add an else statement here like above but with `kc_total_hits = 0`
                bgc_num += 1
            else:
                pass
    # prepares a pandas dataframe from all the dicts
    array = pd.DataFrame(bgc_list)
    if "kc_total_hits" in array:
        array["kc_total_hits"] = array["kc_total_hits"].fillna(0).astype(int)
    
    # code to remove duplicate BGCs resulting from 'neighbouring' BGCs
    # makes a list of all 'neighbouring' clusters
    try:
        grouped = array.groupby('contig')
        to_drop = []
        for group in grouped:
            single_locations = group[group['kind'] == 'single']['BGC_location'].tolist()
            neighbouring_locations = group[group['kind'] == 'neighbouring']['BGC_location'].tolist()
            for single in single_locations:
                single_left, single_right = [int(x) for x in single.split(':')]
                for neighbouring in neighbouring_locations:
                    neighbouring_left, neighbouring_right = [int(x) for x in neighbouring.split(':')]
                    if single_left >= neighbouring_left and single_right <= neighbouring_right:
                        to_drop.extend(group[group['BGC_location'] == single].index)
        array.drop(to_drop, inplace=True)
    except KeyError:
        pass
    
    if len(array) > 0:
        # removes list syntax from 'BGC_type'
        array['BGC_type'] = array['BGC_type'].astype(str).str.replace("[\]\[]",'').str.replace("'",'')
        # assigns an ID to each BGC
        array.insert(1, 'ID', range(1, 1 + len(array)))

    return array

def parse_antismash(input_args, filenames):
    """
    Will open antiSMASH JSON and parse to retrieve information about the number 
    and type of each BGC identified by antiSMASH for a given output folder.
    
    Input:      antismash output folder
    Output:     pandas array with antiSMASH BGC information    
    """

    contents = open_json(filenames.antismash_json)
    df = parse_json(contents, filenames.antismash_json)
    try:
        df.to_csv(filenames.bgc_results)

        flat = {}
        flat["ID"] = input_args.array
        flat["BGC_count"] = str(df.shape[0])
        flat["contig_edge_proportion"] = round(int(df['partial'].value_counts().loc["True"]) / int(df['partial'].size), 2)
        flat["kc_proportion"] = round(len(df[df['kc_mibig_acc'] != '']) / int(df['kc_mibig_acc'].size), 2)
        bgc_list = []
        bgc_list.append(df['BGC_type'].unique().tolist())
        flat["BGC_types"] = bgc_list
        flat["BGC_mean_length"] = round(df['BGC_length'].mean(), 2)
        flat["BGC_median_length"] = round(df['BGC_length'].median(), 2)

        flat_df = pd.DataFrame.from_dict(flat)
        
        return flat_df

    except KeyError:
        print("No BGCs.")

def main(input_args,filenames):

    lib.print_h("Initializing \'post analysis\' module...")
    post_start_time = datetime.datetime.now()
    final_df = pd.DataFrame()
    quast_text = "  _______\n___/ Quast \____________________________________________________________________"
    lib.print_t(quast_text)
    
    results_path = os.path.join(input_args.directory_new,"results")
    quast_path = os.path.join(input_args.directory_new,"quast")
    filenames.quast_report = os.path.join(quast_path,"transposed_report.tsv")
    lib.make_path(quast_path)
    lib.make_path(results_path)
    if lib.file_exists(filenames.quast_report,"Quast already analysed the assembly!","") is False:
        quast(input_args,filenames,quast_path)
    final_df = parse_quast(filenames,input_args.array)

    if input_args.its_db is not None:
        tax_lookup_text = "  ___________________\n___/  Taxonomy Lookup  \________________________________________________________"        
        lib.print_t(tax_lookup_text)

        itsx_path = os.path.join(input_args.directory_new,"ITSx")
        lib.make_path(itsx_path)
        filenames.blast_out = os.path.join(itsx_path,f"{input_args.array}_its_blastn.out")
        its_full = os.path.join(itsx_path,f"{input_args.array}.full.fasta")
        its_2 = os.path.join(itsx_path,f"{input_args.array}.ITS2.fasta")
        its_1 = os.path.join(itsx_path,f"{input_args.array}.ITS1.fasta")
        its_chi = os.path.join(itsx_path,f"{input_args.array}.chimeric.fasta")
        
        if lib.file_exists(its_full,"","") is False and \
            lib.file_exists(its_2,"","") is False and \
                lib.file_exists(its_1,"","") is False and \
                    lib.file_exists(its_chi,"","") is False:
            itsx(input_args,filenames,itsx_path)
        if lib.file_exists(its_full,"ITSx successfully extracted the full ITS sequence!","ITSx failed to extract the full ITS sequence. Checking for ITS2"):
            filenames.its_fasta = its_full
        elif lib.file_exists(its_2,"ITSx successfully extracted the ITS2 sequence!","ITSx failed to extract the ITS2 sequence. Checking for ITS1"):
            filenames.its_fasta = its_2
        elif lib.file_exists(its_1,"ITSx successfully extracted the ITS1 sequence!","ITSx failed to extract the ITS1 sequence. Checking for a chimeric ITS sequence"):
            filenames.its_fasta = its_1
        elif lib.file_exists(its_chi,"ITSx successfully extracted a chimeric ITS sequence!","ITSx failed... check the logs and your inputs"):
            filenames.its_fasta = its_chi
        else:
            lib.print_h("Skipping taxonomy lookup of isolate as there is no extracted ITS sequence!")

        try:    
            if lib.file_exists(filenames.its_fasta,"","") is True:
                blastn(input_args,filenames,itsx_path)
                if lib.file_exists(filenames.blast_out,f"BLASTn successfully searched {filenames.its_fasta} against the ITS_RefSeq database!",f"BLASTn failed to search {filenames.its_fasta}... check the logs and your inputs") is True:
                    its_df = parse_blast(filenames.blast_out,input_args.array)
                    print(its_df.head(1))
                    final_df = final_df.merge(its_df, on="ID")
        except AttributeError:
            pass
    else:
        lib.print_h("Skipping taxonomy lookup of isolate")

    # Run antiSMASH if present in the input arguments
    if input_args.antismash is not None:

        bgc_text = "  ____________\n___/ BGC search \_______________________________________________________________"
        lib.print_t(bgc_text)

        antismash_path = os.path.join(input_args.directory_new,"antismash")
        filenames.antismash_json = os.path.join(antismash_path,f"{input_args.array}.json")
        filenames.antismash_svg = os.path.join(results_path,f"{input_args.array}_bgctype_vs_bgclength.svg")
        filenames.bgc_results = os.path.join(results_path,f"{input_args.array}_bgc.csv")

        # If antiSMASH JSON does not exist, run antiSMASH
        if lib.file_exists(filenames.antismash_json,f"antiSMASH already analysed!",f"Need to run antiSMASH") is False:
            # need to remove existing failed antiSMASH directory otherwise antiSMASH will exit
            # then remake it
            if os.path.exists(antismash_path):
                lib.print_n("Removing existing antiSMASH output directory")
                shutil.rmtree(antismash_path)
            lib.make_path(antismash_path)
            # Will preferentially use the functionally annotated GBK file
            if lib.file_exists(filenames.funannotate_func_gbk,f"Using {filenames.funannotate_func_gbk} as input for antiSMASH", \
                f"Using {filenames.assembly_fasta} as input for antiSMASH") == True:
                filenames.antismash_assembly = filenames.funannotate_func_gbk
            # else will use the prediction GBK file
            elif lib.file_exists(filenames.funannotate_gbk,f"Using {filenames.funannotate_gbk} as input for antiSMASH", \
                f"Using {filenames.assembly_fasta} as input for antiSMASH") == True:
                filenames.antismash_assembly = filenames.funannotate_gbk
            # else will use the sorted assembly with nicer contig names if Funannotate has been run
            elif lib.file_exists(filenames.funannotate_sort_fasta,f"Using {filenames.funannotate_sort_fasta} as input for antiSMASH", \
                f"Using {filenames.assembly_fasta} as input for antiSMASH") == True:
                filenames.antismash_assembly = filenames.funannotate_sort_fasta
            else:
                filenames.antismash_assembly = filenames.assembly_fasta
            antismash(input_args,filenames,antismash_path)
        
        # If the antiSMASH JSON file exists, parse it and plot some figures
        # Might need to add in an exception to catch any errors arising from trying to parse an incomplete antiSMASH run
        if lib.file_exists(filenames.antismash_json,"","") is True:
            bgc_df = parse_antismash(input_args, filenames)
            print(bgc_df)
            try:
                final_df = final_df.merge(bgc_df, on="ID")
                final_df.drop(columns=["ID"])
                plot_strip(bgc_df, input_args, results_path)
                plot_violin(bgc_df, input_args, results_path)
            except TypeError:
                print("No BGCs identified. Skipping BGC plots...")
            except KeyError:
                pass
    else:
        lib.print_n("Skipping BGC search with antiSMASH. Use '--antismash' as a script argument if you would like to do this.")

    lib.print_n("Generating final dataframe")
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):
        lib.print_n(final_df)
    
    filenames.results_csv = os.path.join(results_path,f"{input_args.array}_results.csv")
    final_df.to_csv(filenames.results_csv, index=False)
    lib.file_exists_exit(filenames.results_csv,"Final results report successfully generated","Final results report was not generated")
    lib.print_h(f"post analysis module completed in {datetime.datetime.now() - post_start_time}")

if __name__ == '__main__':
    main()
