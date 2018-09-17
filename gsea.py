# import pdb; pdb.set_trace()
import argparse
from os import listdir, makedirs
from os.path import join, dirname, exists
import urllib.request, urllib.error, urllib.parse
import json

import pandas as pd
from gprofiler import GProfiler
import mygene


REST_URL = "http://data.bioontology.org"
API_KEY = "dc459e3b-6630-4254-949b-20fa522de0a5"

# input: list of genes
# output: list of enriched GO-terms, sorted by p-value


# input: folder containing DE gene lists
# output: in a new folder (input_folder_name_GSEA), gsea results of each list in input

def create_parser():
    parser = argparse.ArgumentParser(description="Gene Set Enrichment Analysis using GO for scQuery's DE gene lists", fromfile_prefix_chars="@")
    parser.add_argument("input_folder", help="Folder which contains the input DE gene lists")
    parser.add_argument("output_folder", help="Output folder to place results in")
    parser.add_argument("-n", "--n_genes", help="Use the top n most significant genes", type=int, default=50)
    parser.add_argument("-g", "--go_branch", help="Which branch of GO to use", choices=["GO:BP", "GO:MF", "GO:CC"], default="GO:BP")
    parser.add_argument("-c", "--correction", help="Algorithm used for multiple testing correction", choices=["GSCS", "FDR", "Bonferroni"], default="FDR")
    parser.add_argument("-r", "--rows", help="Maximum number of ontology terms to print in latex tables", type=int, default=10)
    parser.add_argument("-b", "--background_gene_set", help="Gene set to use as statistical background")
    parser.add_argument("-o", "--ordered", help="Input gene list is ordered (most important first)", action="store_true")
    return parser

def get_correction_method(user_selection):
    ["GSCS", "FDR", "Bonferroni"]
    if user_selection == "GSCS":
        return GProfiler.THR_GSCS
    elif user_selection == "FDR":
        return GProfiler.THR_FDR
    elif user_selection == "Bonferroni":
        return GProfiler.THR_BONFERRONI

def get_json(url):
    opener = urllib.request.build_opener()
    opener.addheaders = [('Authorization', 'apikey token=' + API_KEY)]
    return json.loads(opener.open(url).read())

def get_ontology_term(ontology_ID):
    result = get_json(REST_URL + "/search?require_exact_match=true&ontologies=UBERON,CL&q=" + ontology_ID)["collection"]
    return result[0]['prefLabel']

def write_to_table(cell_type, cell_type_info, gprof_results, filepath, nrows):
    table_tex = ""
    table_tex += r"\begin{longtable}{|l|l|r|}"+"\n"
    table_tex += r"\hline"+"\n"
    # Write header rows
    table_tex += r"\textbf{Cell Type} & \multicolumn{2}{l|}{"+cell_type_info['ontology_term']+r"}\\"+"\n"
    table_tex += r"\hline"+"\n"
    table_tex += r"\textbf{Term ID} & \multicolumn{2}{l|}{"+cell_type+r"}\\"+"\n"
    table_tex += r"\hline"+"\n"
    table_tex += r"\# of experiments & \multicolumn{2}{l|}{"+str(cell_type_info['experiment_count'])
    if cell_type_info['is_pooled']:
        table_tex += " (pooled)"
    table_tex += r"}\\"+"\n"
    table_tex += r"\hhline{|=|=|=|}"+"\n"
    table_tex += r"\textbf{GO ID} & \textbf{Name} & \textbf{p-value} \\" + "\n"
    for i in range(min(nrows, len(gprof_results))):
        table_tex += r"\hline"+"\n"
        term_name = gprof_results[i][11]
        term_name = (term_name[:45] + '...') if len(term_name) > 45 else term_name
        table_tex += r"{} & {} & {:.2e} \\".format(gprof_results[i][8], term_name, float(gprof_results[i][2]))+"\n"
    table_tex += r"\hline"+"\n"
    table_tex += r"\end{longtable}"+"\n"
    with open(filepath, 'w') as f:
        f.write(r"\documentclass{article}"+"\n")
        f.write(r"\usepackage{hhline}"+"\n")
        f.write(r"\usepackage{longtable}"+"\n")
        f.write(r"\usepackage{geometry}"+"\n")
        f.write(r"\geometry{letterpaper, margin=1in}"+"\n")
        f.write(r"\begin{document}"+"\n")
        f.write(table_tex)
        f.write(r"\end{document}"+"\n")
    return table_tex
    
if __name__ == "__main__":
    print("hello...")
    args = create_parser().parse_args()
    if not exists(args.output_folder):
        makedirs(args.output_folder)

    gp = GProfiler("MyTool/0.1", want_header=False)
    
    # Get list of cell types in input, and how many experiments done for each
    cell_types_list = []
    for filename in listdir(args.input_folder):
        if filename.endswith("_meta.csv"):
            cell_types_list.append(filename.split('_')[0])
    cell_type_info = {}
    for cell_type in cell_types_list:
        cell_type_info[cell_type] = {
            'ontology_term': get_ontology_term(cell_type),
            'experiment_count': 0,
            'is_pooled': False,
            'meta_file_path': ""
        }
        for filename in listdir(args.input_folder):
            if filename.startswith(cell_type):
                if filename.endswith("_meta.csv"):
                    cell_type_info[cell_type]['meta_file_path'] = join(args.input_folder, filename)
                elif filename.endswith("_combined.csv"):
                    cell_type_info[cell_type]['is_pooled'] = True
                    cell_type_info[cell_type]['experiment_count'] += 1
                else:
                    cell_type_info[cell_type]['experiment_count'] += 1

    # Get statistical background, if provided
    background_gene_set = None
    if args.background_gene_set is not None:
        print("Using user provided background gene set")
        background_gene_set = pd.Series.from_csv(args.background_gene_set, header=None, index_col=False).tolist()
    
    overall_table_all_f = open(join(args.output_folder, "all_results_all.tex"), 'w')
    overall_table_all_f.write(r"\documentclass{article}"+"\n")
    overall_table_all_f.write(r"\usepackage{hhline}"+"\n")
    overall_table_all_f.write(r"\usepackage{longtable}"+"\n")
    overall_table_all_f.write(r"\usepackage{geometry}"+"\n")
    overall_table_all_f.write(r"\geometry{letterpaper, margin=1in}"+"\n")
    overall_table_all_f.write(r"\begin{document}"+"\n")
    overall_table_up_f = open(join(args.output_folder, "all_results_up.tex"), 'w')
    overall_table_up_f.write(r"\documentclass{article}"+"\n")
    overall_table_up_f.write(r"\usepackage{hhline}"+"\n")
    overall_table_up_f.write(r"\usepackage{longtable}"+"\n")
    overall_table_up_f.write(r"\usepackage{geometry}"+"\n")
    overall_table_up_f.write(r"\geometry{letterpaper, margin=1in}"+"\n")
    overall_table_up_f.write(r"\begin{document}"+"\n")
    overall_table_dn_f = open(join(args.output_folder, "all_results_dn.tex"), 'w')
    overall_table_dn_f.write(r"\documentclass{article}"+"\n")
    overall_table_dn_f.write(r"\usepackage{hhline}"+"\n")
    overall_table_dn_f.write(r"\usepackage{longtable}"+"\n")
    overall_table_dn_f.write(r"\usepackage{geometry}"+"\n")
    overall_table_dn_f.write(r"\geometry{letterpaper, margin=1in}"+"\n")
    overall_table_dn_f.write(r"\begin{document}"+"\n")

    for cell_type, info_dict in cell_type_info.items():
        print("Doing GSEA for {}:".format(cell_type))
        df = pd.read_csv(info_dict['meta_file_path'])
        df_all = df.copy()
        df_all.loc[:,'Avg_log2_fold_change'] = df_all.loc[:,'Avg_log2_fold_change'].abs()
        df_all = df_all.sort_values(by=['Max_adj_p_value', 'Avg_log2_fold_change'], ascending=[True, False])
        #df.to_csv(join(out_path, node_name + "_sorted.csv"))
        print(df_all.shape)
        if df_all.shape[0] < args.n_genes:
            # raise Exception("cell type '{}' has less than {} significant genes!".format(cell_type, args.n_genes))
            print("cell type '{}' has less than {} significant genes!".format(cell_type, args.n_genes))
            continue
        df_all = df_all.head(args.n_genes)
        gene_list_all = df_all['EntrezID'].astype('str').values
        mg = mygene.MyGeneInfo()
        result = mg.getgenes(gene_list_all, fields='symbol', species='mouse')
        gene_list_all = [d['symbol'] for d in result]
        results_all = gp.gprofile(gene_list_all, organism="mmusculus", ordered=args.ordered, correction_method=get_correction_method(args.correction), src_filter=[args.go_branch], custom_bg=background_gene_set)
        print("\t # results returned (all genes) = {}".format(len(results_all)))
        # if len(results) == 0:
        #     continue
        filepath = join(args.output_folder, "{}_all.tex".format(cell_type.replace(':', '_')))
        table_tex = write_to_table(cell_type, info_dict, results_all, filepath, args.rows)
        overall_table_all_f.write(table_tex)
        overall_table_all_f.write("\n")
        filepath = join(args.output_folder, "{}_all_full.tex".format(cell_type.replace(':', '_')))
        write_to_table(cell_type, info_dict, results_all, filepath, len(results_all))

        # Split into Up/Down regulated genes
        df_up = df.loc[df['Avg_log2_fold_change'] > 0]
        df_dn = df.loc[df['Avg_log2_fold_change'] < 0]
        df_dn = df_dn.copy() # Hack to get rid of a warning
        df_dn.loc[:,'Avg_log2_fold_change'] = df_dn.loc[:,'Avg_log2_fold_change'].abs()
        df_up = df_up.sort_values(by=['Max_adj_p_value', 'Avg_log2_fold_change'], ascending=[True, False])
        df_dn = df_dn.sort_values(by=['Max_adj_p_value', 'Avg_log2_fold_change'], ascending=[True, False])
        print(df_up.shape)
        print(df_dn.shape)
        df_up = df_up.head(args.n_genes)
        df_dn = df_dn.head(args.n_genes)
        gene_list_up = df_up['EntrezID'].astype('str').values
        mg = mygene.MyGeneInfo()
        result = mg.getgenes(gene_list_up, fields='symbol', species='mouse')
        gene_list_up = [d['symbol'] for d in result]
        results_up = gp.gprofile(gene_list_up, organism="mmusculus", ordered=args.ordered, correction_method=get_correction_method(args.correction), src_filter=[args.go_branch], custom_bg=background_gene_set)
        print("\t # results returned (up genes) = {}".format(len(results_up)))
        filepath = join(args.output_folder, "{}_up.tex".format(cell_type.replace(':', '_')))
        table_tex = write_to_table(cell_type, info_dict, results_up, filepath, args.rows)
        overall_table_up_f.write(table_tex)
        overall_table_up_f.write("\n")
        filepath = join(args.output_folder, "{}_up_full.tex".format(cell_type.replace(':', '_')))
        write_to_table(cell_type, info_dict, results_up, filepath, len(results_up))
        gene_list_dn = df_dn['EntrezID'].astype('str').values
        mg = mygene.MyGeneInfo()
        result = mg.getgenes(gene_list_dn, fields='symbol', species='mouse')
        gene_list_dn = [d['symbol'] for d in result]
        results_dn = gp.gprofile(gene_list_dn, organism="mmusculus", ordered=args.ordered, correction_method=get_correction_method(args.correction), src_filter=[args.go_branch], custom_bg=background_gene_set)
        print("\t # results returned (dn genes) = {}".format(len(results_dn)))
        filepath = join(args.output_folder, "{}_dn.tex".format(cell_type.replace(':', '_')))
        table_tex = write_to_table(cell_type, info_dict, results_dn, filepath, args.rows)
        overall_table_dn_f.write(table_tex)
        overall_table_dn_f.write("\n")
        filepath = join(args.output_folder, "{}_dn_full.tex".format(cell_type.replace(':', '_')))
        write_to_table(cell_type, info_dict, results_dn, filepath, len(results_dn))

        
    overall_table_all_f.write(r"\end{document}"+"\n")
    overall_table_all_f.close()
    overall_table_up_f.write(r"\end{document}"+"\n")
    overall_table_up_f.close()
    overall_table_dn_f.write(r"\end{document}"+"\n")
    overall_table_dn_f.close()
