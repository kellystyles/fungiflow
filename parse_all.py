import os
import json
import pandas as pd
import library as lib

"""
Standalone script that will parse all outputs and prepare a CSV.
"""

def parse_quast(quast_report):
    """
    Parses `quast` output reports, returning a pandas dataframe.

    Input:      quast TSV report
    Output:     pandas DataFrame object
    """

    lib.print_n(f"Parsing Quast results for {quast_report}")
    try:
        dataframe = pd.read_csv(quast_report, sep="\t", header=0)
        return dataframe
    except pd.errors.EmptyDataError:
        lib.print_e(f"The QUAST report {quast_report} was not parsed successfully")

    return dataframe

def parse_busco(file, final_df):

    with open(file, "r") as f:
        lines = f.readlines()
        
    busco_searched = lines[14].split()[0] 
    print(busco_searched)
    busco_complete = f"{lines[9].split()[0]} ({round(int(lines[9].split()[0])/int(busco_searched) * 100, 2)}%)"
    busco_sc = f"{lines[10].split()[0]} ({round(int(lines[10].split()[0])/int(busco_searched) * 100, 2)}%)"
    busco_dup = f"{lines[11].split()[0]} ({round(int(lines[11].split()[0])/int(busco_searched) * 100, 2)}%)"
    busco_frag = f"{lines[12].split()[0]} ({round(int(lines[12].split()[0])/int(busco_searched) * 100, 2)}%)"
    busco_miss = f"{lines[13].split()[0]} ({round(int(lines[13].split()[0])/int(busco_searched) * 100, 2)}%)"
    print(busco_complete)
    final_df["BUSCO - C"] = [busco_complete] 
    final_df["BUSCO - S"] = [busco_sc] 
    final_df["BUSCO - D"] = [busco_dup]
    final_df["BUSCO - F"] = [busco_frag]
    final_df["BUSCO - M"] = [busco_miss]    

def parse_funannotate(file, final_df):
    contents = open_json(file)
    final_df["genes"] = [contents['annotation']['genes']]
    final_df["mRNA"] = [contents['annotation']['mRNA']]
    final_df["tRNA"] = [contents['annotation']['tRNA']]
    final_df["ncRNA"] = [contents['annotation']['ncRNA']]
    final_df["rRNA"] = [contents['annotation']['rRNA']]
    final_df["avg_gene_length"] = [contents['annotation']['avg_gene_length']]

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

def parse_antismash(json_object,json_file):
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
        x = array.index[array['kind'] == "neighbouring"].tolist()
        for i in x:
            n = x[0]
            for j in range(len(array.loc[i]["BGC_type"])):
                n += 1
                og_left, og_right = array.loc[i]["BGC_location"].split(":")[0], array.loc[i]["BGC_location"].split(":")[1]
                # checks if the BGCs are in the same contig
                if array.loc[i]["contig"] == array.loc[n]["contig"]:
                    # checks if the 'BGC_type' is in the 'neihgbouring' 'BGC_type' list
                    if str(array.loc[n]["BGC_type"][0]) in array.loc[i]["BGC_type"]:
                        # checks if the location is within the range of the original 'location'
                        left, right = array.loc[n]["BGC_location"].split(":")[0], array.loc[n]["BGC_location"].split(":")[1]
                        if left >= og_left and right <= og_right:
                            array.drop([n], axis=0, inplace=True)
    except KeyError:
        pass

def flatten_antismash(antismash_json, bgc_results, final_df):
    """
    Will open antiSMASH JSON and parse to retrieve information about the number 
    and type of each BGC identified by antiSMASH for a given output folder.
    
    Input:      antismash output folder
    Output:     pandas array with antiSMASH BGC information    
    """

    contents = open_json(antismash_json)
    bgc_df = parse_json(contents, antismash_json)
    try:
        if lib.file_exists(bgc_results, "", "") is False:
            bgc_df.to_csv(bgc_results)

        bgc_list = []
        bgc_list.append(bgc_df['BGC_type'].unique().tolist())
        final_df["BGC_count"] = [str(df.shape[0])]
        final_df["contig_edge_proportion"] = [round(int(df['partial'].value_counts().loc["True"]) / int(df['partial'].size), 2)]
        final_df["kc_proportion"] = [round(len(df[df['kc_mibig_acc'] != '']) / int(df['kc_mibig_acc'].size), 2)]
        final_df["BGC_types"] = [bgc_list]
        final_df["BGC_mean_length"] = [round(df['BGC_length'].mean(), 2)]
        final_df["BGC_median_length"] = [round(df['BGC_length'].median(), 2)]

    except KeyError:
        print("No BGCs.")


def main():
    
    master = pd.DataFrame()

    for dir in os.listdir("."):
        os.chdir(dir)
        if lib.file_exists(filenames.funannotate_func_gbk,"Obtaining Funannotate statistics",""):
            if lib.file_exists(os.path.join("funannotate","annotate_misc","run_busco","short_summary_busco.txt"),"",""):
                busco_stats = os.path.join("funannotate","annotate_misc","run_busco","short_summary_busco.txt")
            elif lib.file_exists(os.path.join("funannotate","predict_misc","busco_proteins",f"run_{dir}",f"short_summary_{dir}.txt"),"",""):
                busco_stats = os.path.join("funannotate","annotate_misc","run_busco","short_summary_busco.txt")

        quast_report = os.path.join("quast", "transposed_report.tsv")
        funannotate_json = os.path.join("funannotate", "predict_results", f"{dir}.stats.json")
        antismash_json = os.path.join("antismash", f"{dir}.json")
        bgc_results = os.path.join("results", f"{dir}_bgc.csv")
        outfile = os.path.join("results", f"{dir}_results.csv")
        df = parse_quast(quast_report)
        parse_busco(busco_stats, df)
        parse_funannotate(funannotate_json, df)
        flatten_antismash(antismash_json, bgc_results, df)
        df.to_csv(outfile)
        master = pd.concat([master, df], axis=0)

if __name__ == '__main__':
    main()