'''
Created on 28 Sep 2017

@author: husensofteng
'''

import os, json, sys
from glob import glob

"""Combine data from .bed files listed under data_tracks in the main conf file
if the data was not already combined"""


def collect_all_data(data_dir , data_tracks, sep='\t'):
    
    """Combine all data tracks into a bed4 files one per chr, also record assay types"""
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
        print("Generating chromatin data for all the cells")
        for data_track in data_tracks.split(','):
            print(data_track)
            if os.path.exists(data_track) or '*' in data_track:
                #on linux use awk to generate a file per chr
                if sys.platform == "linux" or sys.platform == "linux2":
                    print("""awk '{print $0 >> \"""" + data_dir + """/"$1".bed"}' """ + data_track)
                    os.system("""awk '{print $0 >> \"""" + data_dir + """/"$1".bed"}' """ + data_track)
                else: #otherwise use readline (it is much slower than awk)
                    for f in glob(data_track):
                        #read the file content and add each line to the chrX file based on col1
                        with open(f, 'r') as fi:
                            l = fi.readline()
                            sl = l.strip().split(sep)
                            while len(sl)>1:
                                with open(data_dir+'/'+sl[0]+'.bed', 'a') as fo:
                                    fo.write(l)
                                l = fi.readline()
                                sl = l.strip().split(sep)
                    
        print("Combined data from the listed tracks.")
    else:
        print("Using existing data tracks from: " + data_dir)
    return data_dir

def get_assay_cell_info(data_dir , 
                        sep, 
                        matching_rep_cell_names_dict, 
                        generated_dicts_output_file, 
                        tissues_with_gene_expression):
    '''Get information (dicts) about the cells and assays from the input data_dir files'''
    assay_cells = {} # e.g DNase-seq: [HepG2,K562,...],...
    assay_cells_datatypes = {}
    cell_assays = {} # e.g HepG2: [DNase-seq,TFBinding,ReplicationDomains], K562: [...], ...
    cell_tfs = {} # e.g HepG2: [CTCF, FOXA1,...], K562: [...], ...
    tf_cells = {} # e.g CTCF: [HepG2, K562,...], FOXA1: [...], ...
    if not os.path.exists(generated_dicts_output_file):
        #add cells/tissues that have gene expr to the dicts
        assay_name = 'TFExpr'
        assay_cells[assay_name] = []
        assay_cells_datatypes[assay_name] = 'real'
        for tissue in tissues_with_gene_expression:
            assay_cells[assay_name].append(tissue)
            cell_assays[tissue] = [assay_name]
        
        for chr_file in os.listdir(data_dir):
            if "chr" in chr_file:
                print('reading from: ' + data_dir + '/' + chr_file)
                with open(data_dir+'/'+chr_file, 'r') as data_infile:
                    l = data_infile.readline()
                    while l:
                        sl = l.strip().split(sep)
                        if len(sl)<4:#there must be at least 4 cols
                            l = data_infile.readline()
                            continue
                        tracks = sl[3].split(',')
                        for track in tracks:
                            info_col = track.split('#')
                            cell_name = info_col[0]
                            rep_cell_names = [cell_name]
                            assay = info_col[1]
                            #get the datatype it can be used for databse storage
                            if len(info_col)==2:
                                assay_cells_datatypes[assay] = 'real'
                            elif info_col[1]=="TFBinding":
                                assay_cells_datatypes[assay] = 'real'
                                assay_cells_datatypes['NumOtherTFBinding'] = 'real'
                                assay_cells_datatypes['OtherTFBinding'] = 'text'
                            else:
                                try:
                                    float(info_col[-1])
                                    assay_cells_datatypes[assay] = 'real'
                                except ValueError:
                                    assay_cells_datatypes[assay] = 'text'
                            
                            try:#get the cell's rep-name
                                rep_cell_names = matching_rep_cell_names_dict[cell_name]
                            except KeyError:
                                rep_cell_names = [info_col[0]]
                            for cell in rep_cell_names:
                                try:#add the cell name to the corresponding assay list in assay_cells
                                    if cell not in assay_cells[assay]: 
                                        assay_cells[assay].append(cell)
                                except KeyError:
                                    assay_cells[assay] = [cell]
                                if assay=="TFBinding":
                                    try:#add the cell name to the corresponding assay list in assay_cells
                                        if cell not in assay_cells["NumOtherTFBinding"]: 
                                            assay_cells["NumOtherTFBinding"].append(cell)
                                    except KeyError:
                                        assay_cells["NumOtherTFBinding"] = [cell]
                                    try:
                                        if cell not in assay_cells["OtherTFBinding"]: 
                                            assay_cells["OtherTFBinding"].append(cell)
                                    except KeyError:
                                        assay_cells["OtherTFBinding"] = [cell]
                            
                            for cell in rep_cell_names:
                                try:#add the assay name to the corresponding cell list in cells_assay
                                    if assay not in cell_assays[cell]: 
                                        cell_assays[cell].append(assay)
                                except KeyError:
                                    cell_assays[cell] = [assay]
                                if assay=="TFBinding":    
                                    try:
                                        if "NumOtherTFBinding" not in cell_assays[cell]: 
                                            cell_assays[cell].append("NumOtherTFBinding")
                                    except KeyError:
                                        cell_assays[cell] = ["NumOtherTFBinding"]
                                    try:
                                        if "OtherTFBinding" not in cell_assays[cell]: 
                                            cell_assays[cell].append("OtherTFBinding")
                                    except KeyError:
                                        cell_assays[cell] = ["OtherTFBinding"]
                                        
                                    
                            if assay=="TFBinding":
                                tf_name = info_col[2]
                                for cell in rep_cell_names:
                                    try:
                                        if tf_name not in cell_tfs[cell]: 
                                            cell_tfs[cell].append(tf_name)
                                    except KeyError:
                                        cell_tfs[cell] = [tf_name]
                                    try:
                                        if cell not in tf_cells[tf_name]: 
                                            tf_cells[tf_name].append(cell)
                                    except KeyError:
                                        tf_cells[tf_name] = [cell]
                        l = data_infile.readline()
        
        #write the dict contents to an output file
        with open(generated_dicts_output_file, 'w') as generated_dicts_outfile:
            generated_dicts_outfile.write("assay_cells=")
            json.dump(assay_cells, generated_dicts_outfile)
            generated_dicts_outfile.write('\ncell_assays=')
            json.dump(cell_assays, generated_dicts_outfile)
            generated_dicts_outfile.write('\ncell_tfs=')
            json.dump(cell_tfs, generated_dicts_outfile)
            generated_dicts_outfile.write('\ntf_cells=')
            json.dump(tf_cells, generated_dicts_outfile)
            generated_dicts_outfile.write("\nassay_cells_datatypes=")
            json.dump(assay_cells_datatypes, generated_dicts_outfile)
    else:
        with open(generated_dicts_output_file, 'r') as generated_dicts_outfile:
            lines_from_dict_file = generated_dicts_outfile.readlines()
            number_of_loaded_dicts = 0
            for l in lines_from_dict_file:
                if '=' in l:
                    if l.split('=')[0]=="assay_cells":
                        assay_cells = json.loads(l.strip().split('=')[1])
                        number_of_loaded_dicts+=1
                    elif l.split('=')[0]=="cell_assays":
                        cell_assays = json.loads(l.strip().split('=')[1])
                        number_of_loaded_dicts+=1
                    elif l.split('=')[0]=="cell_tfs":
                        cell_tfs = json.loads(l.strip().split('=')[1])
                        number_of_loaded_dicts+=1
                    elif l.split('=')[0]=="tf_cells":
                        tf_cells = json.loads(l.strip().split('=')[1])
                        number_of_loaded_dicts+=1
                    elif l.split('=')[0]=='assay_cells_datatypes':
                        assay_cells_datatypes = json.loads(l.strip().split('=')[1])
                        number_of_loaded_dicts+=1
                    else:
                        continue
            if number_of_loaded_dicts<5:
                sys.exit('Error: the given dict file is corrupted please remove it and re-run the program to regenerate it: ref_file: ' + generated_dicts_output_file)
    #write the content of assay_cells, cell_assays, cell_tfs, tf_cells, assay_cells_datatypes to a file with proper headers for each and re-construct of each from the file if it exists
    return assay_cells, cell_assays, cell_tfs, tf_cells, assay_cells_datatypes

def generate_cells_assays_matrix(cell_assays, 
                                 cell_names, 
                                 assay_cells_datatypes, 
                                 tissues_with_gene_expression):
    '''Generate a default dict based on the information obtained from the tracks in data_dir'''
    cells_assays_dict = {}
    for cell_name in cell_assays.keys():
        if not (cell_name in cell_names or cell_name in tissues_with_gene_expression):#only consider cells that are listed in the dict or the gene expr file
            continue
        cells_assays_dict[cell_name] = {}
        for assay_name in cell_assays[cell_name]:
            try:
                if assay_cells_datatypes[assay_name] == "real":
                    cells_assays_dict[cell_name][assay_name] = 0.0
                else:
                    cells_assays_dict[cell_name][assay_name] = "NO"
                if assay_name=="TFBinding":
                    cells_assays_dict[cell_name]["NumOtherTFBinding"] = 0.0
                    cells_assays_dict[cell_name]["OtherTFBinding"] = []
            except ValueError:#consider it as text if no data type was found
                cells_assays_dict[cell_name][assay_name] = "NO"
    return cells_assays_dict

