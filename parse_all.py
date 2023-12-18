import os
import sys
import json
import pandas as pd
import library as lib

"""
Standalone script that will parse all Fungiflow outputs and prepare a CSV.
Input the parent directory containing all the individual Fungiflow directories
as a paramter, otherwise will be assumed to be the current directory. Make sure
that no non-fungiflow directories are in the target directory.

Usage: python3 parse_all.py 'directory'
"""

def parse_quast(quast_report):
    """
    Parses `quast` output reports, returning a pandas dataframe.

    Input:      quast TSV report
    Output:     pandas DataFrame object
    """

    try:
        dataframe = pd.read_csv(quast_report, sep="\t", header=0)
        return dataframe
    except pd.errors.EmptyDataError:
        lib.print_e(f"The QUAST report {quast_report} was not parsed successfully")

def parse_itsx(results_path, final_df):
    """
    Parses BLAST output files and appends the extracted information to a pandas DataFrame.

    Input:      BLAST output file from blast[x]() function [path]
    Output:     pandas DataFrame object
    """

    columns = ["Query","Species","Accession","Identity","Bitscore","E-value","Subject Name"]

    if os.path.exists(results_path):
        try:
            dataframe = pd.read_csv(results_path, sep="\t", header=None)
            dataframe.columns = columns
            dataframe = dataframe.head(1)
            final_df["Query"] = dataframe["Query"]
            final_df["Species"] = dataframe["Species"]
            final_df["Accession"] = dataframe["Accession"]
            final_df["Identity"] = dataframe["Identity"]
            final_df["Bitscore"] = dataframe["Bitscore"]
            final_df["E-value"] = dataframe["E-value"]
            final_df["Subject Name"] = dataframe["Subject Name"]
            return dataframe
        except pd.errors.EmptyDataError:
            print(f"No BLASTp results obtained for {results_path}")

def parse_busco(file, final_df):

    with open(file, "r") as f:
        lines = f.readlines()

    for line in lines:
        if "Total BUSCO groups searched" in line:
            busco_searched = line.split()[0]            
    for line in lines:         
        if "Complete BUSCOs" in line:
            busco_complete = f"{line.split()[0]} ({round(int(line.split()[0])/int(busco_searched) * 100, 2)}%)"
        if "Complete and single-copy BUSCOs" in line:
            busco_sc = f"{line.split()[0]} ({round(int(line.split()[0])/int(busco_searched) * 100, 2)}%)"
        if "Complete and duplicated BUSCOs" in line:
            busco_dup = f"{line.split()[0]} ({round(int(line.split()[0])/int(busco_searched) * 100, 2)}%)"
        if "Fragmented BUSCOs" in line:
            busco_frag = f"{line.split()[0]} ({round(int(line.split()[0])/int(busco_searched) * 100, 2)}%)"
        if "Missing BUSCOs " in line:
            busco_miss = f"{line.split()[0]} ({round(int(line.split()[0])/int(busco_searched) * 100, 2)}%)"

    final_df["BUSCO - C"] = [busco_complete] 
    final_df["BUSCO - S"] = [busco_sc] 
    final_df["BUSCO - D"] = [busco_dup]
    final_df["BUSCO - F"] = [busco_frag]
    final_df["BUSCO - M"] = [busco_miss]

def parse_funannotate(file, final_df):
    contents = open_json(file)
    try:
        final_df["genes"] = [contents['annotation']['genes']]
        final_df["mRNA"] = [contents['annotation']['mRNA']]
        final_df["tRNA"] = [contents['annotation']['tRNA']]
        final_df["ncRNA"] = [contents['annotation']['ncRNA']]
        final_df["rRNA"] = [contents['annotation']['rRNA']]
        final_df["avg_gene_length"] = [contents['annotation']['avg_gene_length']]
    except TypeError:
        pass

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
                return None
    except:
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

def get_chemclass(list):

    n = ["cdps", "nrps", "nrps-like", "thioamide-nrp"]
    p = ["hgle-ks","pks-like","ppys-ks","t1pks","t2pks","t3pks","transat-pks",'transat-pks-like']
    r = ["bottromycin","cyanobactin","fungal-ripp","glycocin","lap","lantipeptide class i","lantipeptide class ii","lantipeptide class iii","lantipeptide class iv","lantipeptide class v","lassopeptide","linaridin","lipolanthine","microviridin","proteusin","ranthipeptide","ras-ripp","ripp-like","rre-containing","sactipeptide","thioamitides","thiopeptide","bacteriocin","head_to_tail","lanthidin","lanthipeptide","tfua-related","microcin","lanthipeptide"]
    s = ["amglyccycl","oligosaccharide","saccharide","cf_saccharide"]
    t = ["terpene"]

    if len(list) == 1:
        if list[0].lower() in n:
            chemclass = "NRP"
        elif list[0].lower() in p:
            chemclass = "Polyketide"
        elif list[0].lower() in r:
            chemclass = "RiPP"
        elif list[0].lower() in s:
            chemclass = "Saccharide"
        elif list[0].lower() in t:
            chemclass = "Terpene"
        else:
            chemclass="Other"            
    if len(list) > 1:
        chemclass="Hybrid"
    return chemclass

def parse_antismash(json_object,json_file):
    """
    Extracts information about the BGCs from a JSON object,
    including knownclusterblast data.
    Input: JSON object (from open_json), JSON file (str)
    Output: pandas array of information
    """
    bgc_list = []
    bgc_num = 1
    input_file = str(json_file.split("/")[-1])
    try:
        for i in range(len(json_object['records'])):
            contig = json_object['records'][i]['id']#.split("_")[0]
            # loops through candidate clusters and extract information about each cluster
            for j in range(len(json_object['records'][i]['features'])):
                record = json_object['records'][i]['features'][j]
                kc_record = json_object['records'][i]["modules"]
                if 'cand_cluster' in record['type']:
                    # prepares dict and adds information to this dict
                    a = {}
                    bgc_id = contig + "_" + str(j)
                    product = record["qualifiers"]['product']
                    edge = record["qualifiers"]["contig_edge"][0]
                    bgc_begin, bgc_end = get_location(record["location"])
                    bgc_length = bgc_end - bgc_begin
                    a["file"] = input_file
                    a["contig"] = contig
                    a["partial"] = edge
                    a["BGC_type"] = product
                    a["chem_class"] = get_chemclass(product)
                    a["BGC_location"] = str(record["location"]).lstrip("[<").rstrip(">]")
                    a["BGC_length"] = bgc_length
                    a["protocluster"] = record["qualifiers"]["protoclusters"][0]
                    a["kind"] = record["qualifiers"]["kind"][0]
                    # extracts the known cluster information
                    if "antismash.modules.clusterblast" in kc_record:
                        for c in kc_record["antismash.modules.clusterblast"]['general']['results']:
                            if len(c['ranking']) != 0:
                                a["cblast_type"] = c['ranking'][0][0]["cluster_type"]
                                a["cblast_accession"] = c['ranking'][0][0]["accession"]
                                a["cblast_tags"] = c['ranking'][0][0]["tags"]
                                genes = len(c['ranking'][0][0]['tags'])
                                temp_l1 = []
                                for item in c['ranking'][0][1]['pairings']:
                                    temp_l1.append(item[2]['name'])
                                cblaster = len(list(set(temp_l1)))*100 / genes
                                if c['ranking'][0][1]["core_gene_hits"] == 0:
                                    item = (bgc_id, "", "")
                                else:
                                    item = (bgc_id, a["cblast_accession"] + ".0", cblaster)
                                #print(a["cblast_type"], a["cblast_accession"], a["cblast_tags"], item)
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
                    bgc_num += 1
                else:
                    pass
    except TypeError:
        pass    
    # prepares a pandas dataframe from all the dicts
    array = pd.DataFrame(bgc_list)
    if "kc_total_hits" in array:
        array["kc_total_hits"] = array["kc_total_hits"].fillna(0).astype(int)
    
    # code to remove duplicate BGCs resulting from 'neighbouring' BGCs
    # makes a list of all 'neighbouring' clusters
    # Group the rows by the 'contig' column

    try:
        grouped = array.groupby('contig')
        to_drop = []
        for group in grouped:
            #print(type(group))
            single_locations = group[1][group[1]['kind'] == 'single']['BGC_location'].tolist()
            neighbouring_locations = group[1][group[1]['kind'] == 'neighbouring']['BGC_location'].tolist()
            for single in single_locations:
                single = single.replace('<', '').replace('>', '')
                single_left, single_right = [int(x) for x in single.split(':')]
                for neighbouring in neighbouring_locations:
                    neighbouring_left, neighbouring_right = [int(x) for x in neighbouring.split(':')]
                    if single_left >= neighbouring_left and single_right <= neighbouring_right:
                        to_drop.extend(group[1][group[1]['BGC_location'] == single].index)
        array.drop(to_drop, inplace=True)
    except KeyError:
        pass

    return array

def flatten_antismash(antismash_json, bgc_results, final_df):
    """
    Will open antiSMASH JSON and parse to retrieve information about the number 
    and type of each BGC identified by antiSMASH for a given output folder.
    
    Input:      antismash output folder
    Output:     pandas array with antiSMASH BGC information    
    """

    contents = open_json(antismash_json)
    bgc_df = parse_antismash(contents, antismash_json)

    try:
        #if lib.file_exists(bgc_results, "", "") is False:
        bgc_df.to_csv(bgc_results)
        flattened_list = [item for sublist in bgc_df['BGC_type'] for item in sublist]
        unique_bgc_class = set(flattened_list)
        final_df["BGC_count"] = [str(bgc_df.shape[0])]
        final_df["kc_proportion"] = [round(len(bgc_df[bgc_df['kc_mibig_acc'] != '']) / int(bgc_df['kc_mibig_acc'].size), 2)]
        final_df["BGC_types"] = ", ".join(unique_bgc_class)
        final_df["total_BGC_length"] = bgc_df['BGC_length'].sum()
        final_df["BGC_mean_length"] = [round(bgc_df['BGC_length'].mean(), 2)]
        final_df["BGC_median_length"] = [round(bgc_df['BGC_length'].median(), 2)]
    except KeyError:
        pass
    try:
        final_df["contig_edge_proportion"] = [round(int(bgc_df['partial'].value_counts().loc["True"]) / int(bgc_df['partial'].size), 2)]
    except KeyError:
        final_df["contig_edge_proportion"] = 0
    
    return bgc_df

def main():
    try:
        main_dir = sys.argv[1]
    except IndexError:
        main_dir = "."

    lib.print_h("|-|-|-|-|-|-|-|- parse_all.py -|-|-|-|-|-|-|-|")
    lib.print_n(f"Parsing all Fungiflow directories in {main_dir}")
    master = pd.DataFrame()
    bgc_master = pd.DataFrame()           
    directories = [d for d in os.listdir(main_dir) if os.path.isdir(os.path.join(main_dir, d))]
    directories[:] = [d for d in directories if os.path.isdir(os.path.join(main_dir, d, "assembly")) is True]
    n = 0

    # Loop through all fungiflow directories and parse information
    for dir in directories:
        n += 1
        lib.print_n(f"Parsing results for {dir:5} - {n:5}/{len(directories)}")
        df = pd.DataFrame()
        try:
            quast_report = os.path.join(main_dir, dir, "quast", "transposed_report.tsv")
            qdf = parse_quast(quast_report)
            df = pd.concat([df, qdf], axis=0)
        except FileNotFoundError:
            # if no Quast files present
            pass  
        try:
            funannotate_json = os.path.join(main_dir, dir, "funannotate", "predict_results", f"{dir}.stats.json")
            parse_funannotate(funannotate_json, df)
        except FileNotFoundError:
            # if no funannotate files present
            pass
        try:
            if os.path.exists(os.path.join(main_dir, dir, "results")):
                bgc_results = os.path.join(main_dir, dir, "results", f"{dir}_bgc.csv")
            else:
                bgc_results = os.path.join(main_dir, dir, f"{dir}_bgc.csv")
            antismash_json = os.path.join(main_dir, dir, "antismash", f"{dir}.json")
            if os.path.exists(antismash_json):
                bgcs = flatten_antismash(antismash_json, bgc_results, df)
                bgc_master = pd.concat([bgc_master, bgcs], axis=0)
        except FileNotFoundError:
            # if no antiSMASH dir present
            pass
        try:
            if lib.file_exists(os.path.join(main_dir, dir, "ITSx",f"{dir}_its_blastn.out"),"",""):
                itsx_output = os.path.join(main_dir, dir, "ITSx",f"{dir}_its_blastn.out")
                parse_itsx(itsx_output, df)
        except FileNotFoundError:
            # if no ITSx BLSTn file present
            pass
        try:

            if lib.file_exists(os.path.join(main_dir, dir, "funannotate","annotate_misc","run_busco","short_summary_busco.txt"),"",""):
                busco_stats = os.path.join(main_dir, dir, "funannotate","annotate_misc","run_busco","short_summary_busco.txt")
                parse_busco(busco_stats, df)
            elif lib.file_exists(os.path.join(main_dir, dir, "funannotate","predict_misc","busco_proteins",f"run_{dir}",f"short_summary_{dir}.txt"),"",""):
                busco_stats = os.path.join(main_dir, dir, "funannotate","predict_misc","busco_proteins",f"run_{dir}",f"short_summary_{dir}.txt")
                parse_busco(busco_stats, df)
        except FileNotFoundError:
            # if no BUSCO files present
            pass

        if os.path.exists(os.path.join(main_dir, dir, "results")):
            outfile = os.path.join(main_dir, dir, "results", f"{dir}_results.csv")
        else:
            outfile = os.path.join(main_dir, dir, f"{dir}_results.csv")
        
        df["Assembly"] = str(dir)
        df.to_csv(outfile)
        master = pd.concat([master, df], axis=0)

    # Saving dataframes to CSVs
    bgc_master.to_csv(os.path.join(main_dir, "all_bgcs.csv"), index=False)
    master.to_csv(os.path.join(main_dir, "master_results.csv"), index=False)
    lib.print_h(f"{len(master)} Fungiflow output directories parsed!")

if __name__ == '__main__':
    main()
