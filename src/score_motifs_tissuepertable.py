'''
Created on Nov 13, 2016

@author: Husen M. Umer

Score motifs: collects cell-type specific data from several public resources and generates a cell-type specific score for each motif instance in the human genome
Input: TF PWMs, human genome, TF chip-seq resources, DNase1 resources, ChromHMM labels, Gene expression, CAGE peaks, HIC domains, HIC loops, Replication domains 
Output: A list of motif instances with a functionality score per cell type 
Process: the module has three sections 1)collects and processes data from the provided resources, 2) combines data from the collections and overlays them with motifs, 3) computes a score for each motif instance
'''
import os
import sys
from collections import Counter
from itertools import islice
import multiprocessing as mp
from multiprocessing import Pool
import json
import psycopg2
import time
import math
from psycopg2.extensions import ISOLATION_LEVEL_AUTOCOMMIT
from psycopg2.extras import DictCursor
from pybedtools import BedTool, set_tempdir
set_tempdir('../src/tmp')

def get_params(params_list):
    params = {}
    for arg in params_list:#priority is for the command line
        if '=' in arg: 
            if len(arg.strip().split('='))==2:
                if arg.split('=')[0] not in list(params.keys()):
                    params[arg.strip().split('=')[0]] = arg.strip().split('=')[1]
    if 'param_file' in params:     
        with open(params['param_file'], 'r') as params_infile:
            params_from_file = params_infile.readlines()
            for line in params_from_file:
                if not line.startswith('//') and not line.startswith('#') and '=' in line:
                    if len(line.strip().split('='))==2:
                        if line.strip().split('=')[0] not in list(params.keys()):
                            params[line.strip().split('=')[0]] = line.strip().split('=')[1].split('#')[0].split('//')[0]
    return params
#change back after test
updated_s_for_test = sys.argv
updated_s_for_test.append('param_file=/data2/husen/ActiveMotifs/score_motifs_param_regulatorymotifs.conf')
params = get_params(updated_s_for_test)

def get_value(str):
    if 'true' in str.lower() or 'yes' in str.lower():
        return True
    else:
        return False

#Section 1: Collect resources
def collect_all_data(data_dir = params['all_chromatin_makrs_all_cells_combined_dir_path'], data_tracks = params['data_tracks']):
    """Combine all data tracks into a bed4 files one per chr, also record assay types"""
    if not os.path.exists(data_dir):
        os.mkdir(data_dir)
        print("Generating chromatin data for all the cells")
        for data_track in data_tracks.split(','):
            print(data_track)
            if os.path.exists(data_track) or '*' in data_track:
                os.system("""awk '{print $0 >> \"""" + data_dir + """/"$1".bed"}' """ + data_track)
        print("Combined data from the listed tracks.")
    else:
        print(("Using existing data tracks from: " + data_dir))
    return data_dir

def get_assay_cell_info(data_dir = params['all_chromatin_makrs_all_cells_combined_dir_path'], sep='\t', matching_rep_cell_names_dict={}, generated_dicts_output_file=params['all_chromatin_makrs_all_cells_combined_dir_path']+"_generated_dicts.txt", tissues_with_gene_expression = {}):
    
    assay_cells = {}#e.g DNase-seq: [HepG2,K562,...],...
    assay_cells_datatypes = {}
    cell_assays = {}#e.g HepG2: [DNase-seq,TFBinding,ReplicationDomains], K562: [...], ...
    cell_tfs = {}#e.g HepG2: [CTCF, FOXA1,...], K562: [...], ...
    tf_cells = {}#e.g CTCF: [HepG2, K562,...], FOXA1: [...], ...
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
                print(('reading from: ' + data_dir + '/' + chr_file))
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

def generate_cells_assays_matrix(cell_assays, cell_names, assay_cells_datatypes, tissues_with_gene_expression):
    cells_assays_dict = {}
    for cell_name in list(cell_assays.keys()):
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

def reset_cells_assays_matrix(tf_name_from_motif_name, cells_assays_dict, cell_tfs, tf_cells, motifTFName_TFNames_matches_dict, assay_cells_datatypes):
    
    for representative_cell in cells_assays_dict:
        for assay in cells_assays_dict[representative_cell]:
            if "TFBinding" not in assay and assay!="TFExpr":
                try:
                    if assay_cells_datatypes[assay] == "real":
                        cells_assays_dict[representative_cell][assay] = 0.0
                    else:
                        cells_assays_dict[representative_cell][assay] = "NO"
                except ValueError:
                    cells_assays_dict[representative_cell][assay] = "NO"
            
            elif "TFExpr" in assay:
                cells_assays_dict[representative_cell][assay] = "NaN"
              
            elif "TFBinding" in assay:#checking wether this TF is available in the current cell
                tf_exists = 0
                if tf_name_from_motif_name in tf_cells:
                    if representative_cell in tf_cells[tf_name_from_motif_name]:
                        cells_assays_dict[representative_cell][assay] = 0.0
                        tf_exists = 1
                else:
                    for alt_tf_name in motifTFName_TFNames_matches_dict[tf_name_from_motif_name]:
                        if alt_tf_name in tf_cells:
                            if representative_cell in tf_cells[alt_tf_name]:
                                cells_assays_dict[representative_cell][assay] = 0.0
                                tf_exists = 1
                                break
                if tf_exists==0:
                    cells_assays_dict[representative_cell][assay]='NaN'
                if len(cell_tfs[representative_cell])-tf_exists > 0.0:
                    cells_assays_dict[representative_cell]["NumOtherTFBinding"] = 0.0
                    cells_assays_dict[representative_cell]["OtherTFBinding"] = []
                else:
                    cells_assays_dict[representative_cell]["NumOtherTFBinding"] = 'NaN'
                    cells_assays_dict[representative_cell]["OtherTFBinding"] = []
    return cells_assays_dict

#Get dicts
def retreive_key_values_from_dict_file(dict_input_file, key_value_sep='=', values_sep=','):#TFFamilyName TF_name
    "Retrieves the key and its values"
    key_values_dict = {}
    value_key_dict = {}
    with open(dict_input_file, 'r') as dict_input_file_infile:
        lines = dict_input_file_infile.readlines()
        for line in lines:
            if line.startswith('//') or line.startswith('#'):# or '=' not in line:
                continue
            sl = line.strip().split(key_value_sep)
            key_value = sl[0].strip()
            if key_value not in list(key_values_dict.keys()):
                key_values_dict[key_value] = []
            if key_value not in value_key_dict:
                value_key_dict[key_value]=[key_value]
            if len(sl)>1:
                for s in sl[1].split(values_sep):
                    if s.strip()!="" and s.strip() not in key_values_dict[key_value]:
                        key_values_dict[key_value].append(s.strip())
                    if s.strip()!="":
                        if s.strip() not in value_key_dict:
                            value_key_dict[s.strip()]=[]
                        value_key_dict[s.strip()].append(key_value)
    return key_values_dict, value_key_dict
    
def retreive_TFFamilyName_for_motifNames():#TFFamilyName TF_name
    "Retrieves the TF family name for each TF name"
    
    motifTFName_TFNames_matches_dict = {}
    with open(params['TF_family_matches_file'], 'r') as TFFamily_matches_infile:
        lines = TFFamily_matches_infile.readlines()
        for line in lines:
            sl = line.strip().split('\t')
            if len(sl)>1:
                motifTF_name = sl[0].strip().upper()
                if motifTF_name not in list(motifTFName_TFNames_matches_dict.keys()):
                    motifTFName_TFNames_matches_dict[motifTF_name] = []
                for s in sl:
                    if s.strip().upper()!="": 
                        if s.strip().upper() not in motifTFName_TFNames_matches_dict[motifTF_name]:
                            motifTFName_TFNames_matches_dict[motifTF_name].append(s.strip().upper())                       
    return motifTFName_TFNames_matches_dict

#Given a GTEX file retrieve gene expression from each tissue for each TF name (as defined in the retreive_TFFamilyName_for_motifNames) 

def get_expression_level_per_originType_per_TF(motifTFName_TFNames_matches_dict, normal_gene_expression_inputfile=params['normal_gene_expression_inputfile'],
                                               origin_gene_expression_values_outputfile = params['normal_gene_expression_inputfile'] + "_perTissue_perTF", 
                                               index_tissues_names_row_start = 2, index_gene_names_col = 1, index_gene_values_start=2, sep='\t'):
    tissue_origin_gene_expression_values = {}
    if not os.path.exists(origin_gene_expression_values_outputfile):
        tf_names_to_extract_gene_expression_for = []#list_tf_names_from_tracks#get names of TFs from the TFFamily file and the dirs contaning ChIP-seq datasets
        for k in list(motifTFName_TFNames_matches_dict.keys()):
            tf_names_to_extract_gene_expression_for.append(k)
        tf_names_to_extract_gene_expression_for = list(set(tf_names_to_extract_gene_expression_for))
        TFs_extract_expression = get_TFs_extract_expression(tf_names_to_extract_gene_expression_for)
        with open(normal_gene_expression_inputfile) as normal_gene_expression_infile:
            for i in range(0, index_tissues_names_row_start):
                line = normal_gene_expression_infile.readline()#skip until it gets to the tissue names row
            tissue_names = normal_gene_expression_infile.readline().strip().split(sep)[index_gene_names_col+1::]
            line = normal_gene_expression_infile.readline()
            
            #iniatilze the tissue_origin_gene_expression_values matrix
            for tissue in tissue_names:
                tissue_origin_gene_expression_values[tissue] = {}
                for tf in TFs_extract_expression:
                    tissue_origin_gene_expression_values[tissue][tf] = 'NaN'
            while line:
                gene_name = line.strip().split(sep)[index_gene_names_col].upper()
                gene_values = line.strip().split(sep)[index_gene_names_col+1::]
                motif_tf_names_corresponding_to_this_gene = get_tfName_fromGeneName(TFs_extract_expression, gene_name)
                if len(motif_tf_names_corresponding_to_this_gene)>0:
                    for tf_name in motif_tf_names_corresponding_to_this_gene:
                        for i in range(0, len(tissue_names)):
                            if tissue_origin_gene_expression_values[tissue][tf_name] == 'NaN':
                                tissue_origin_gene_expression_values[tissue_names[i]][tf_name] = float(gene_values[i])
                            else:
                                tissue_origin_gene_expression_values[tissue_names[i]][tf_name] += float(gene_values[i])
                line = normal_gene_expression_infile.readline()
        #write the results to output file
        with open(origin_gene_expression_values_outputfile, 'w') as origin_gene_expression_values_outfile:
            json.dump(tissue_origin_gene_expression_values,origin_gene_expression_values_outfile)
            #for tissue_name in tissue_origin_gene_expression_values:
            #    for gene_name in tissue_origin_gene_expression_values[tissue_name]:
            #        origin_gene_expression_values_outfile.write(tissue_name + sep + gene_name + sep + str(tissue_origin_gene_expression_values[tissue_name][gene_name]) + '\n')
    else:
        with open(origin_gene_expression_values_outputfile, 'r') as origin_gene_expression_values_infile:
            tissue_origin_gene_expression_values = json.load(origin_gene_expression_values_infile)
    return tissue_origin_gene_expression_values

def get_tfName_fromGeneName(TFs_extract_expression, gene_name):
    motif_tf_names_corresponding_to_this_gene = []
    for tf in TFs_extract_expression:
        if gene_name==tf or gene_name in TFs_extract_expression[tf]:
            motif_tf_names_corresponding_to_this_gene.append(tf)
    return motif_tf_names_corresponding_to_this_gene

def get_TFs_extract_expression(list_of_tf_names):
    TFs_extract_expression = {}
    for tf_name in list_of_tf_names:
        if tf_name not in TFs_extract_expression:
            TFs_extract_expression[tf_name] = []
            if '(' in tf_name:
                TFs_extract_expression[tf_name].append(tf_name.split('(')[0])
            #add TF names after removing dashes and splitting them by :: (in case they were in the name)
            if '-' in tf_name:#in case - was in the name remove it
                TFs_extract_expression[tf_name].append(tf_name.replace('-',''))
            if '::' in tf_name:
                s_tf_name = tf_name.split('::')
                for s in s_tf_name:
                    if s not in TFs_extract_expression[tf_name]:
                        TFs_extract_expression[tf_name].append(s)
    return TFs_extract_expression

#Section2 Overlap between the generated resources and motifs

def run_overlay_resources_score_motifs(normal_expression_per_tissue_origin_per_TF,
                                       matching_cell_name_representative_dict, motifTFName_TFNames_matches_dict, 
                                       cells_assays_dict, cell_tfs, tf_cells, assay_cells_datatypes, header):
    """pairs matching chromosomes in motif_sites_input_dir and all_chromatin_makrs_all_cells_input_dir and calls overlay_resources_score_motifs
    Input: moitf instances input dir (one file per chr)
           chromatin data collection dir (one file per chr, bed4 format; track pos, track cell#assaytype#value or cell#TFname in case of chip-seq) 
    Return: a list of motif_overlapping_track files
    Precondition: files in motif_sites_input_dir and chromatin_tracks_input_dir should have the same names 
                  Recommended: name files in both dirs as chrNumber, chrX or chrY (where number is between 1-22)
    """
    motif_files = []
    if not os.path.isdir(params['motif_sites_dir']) and os.path.isfile(params['motif_sites_dir']):
        motif_files = [params['motif_sites_dir']]
        params['motif_sites_dir'] = "."
    else:
        motif_files = os.listdir(params['motif_sites_dir'])
    
    chromatin_tracks_files = os.listdir(params['all_chromatin_makrs_all_cells_combined_dir_path'])
    if not os.path.exists(params['motifs_overlapping_tracks_output_dir']):
        os.mkdir(params['motifs_overlapping_tracks_output_dir'])
    motifs_overlapping_tracks_files = []
    scored_motifs_overlapping_tracks_files = []
    if get_value(params['run_in_parallel_param']) and len(motif_files)>1:
        p = Pool(int(params['number_processes_to_run_in_parallel']))
    for motif_file in motif_files:
        if motif_file.split('/')[-1] in chromatin_tracks_files:#it is assumed for every motif file name there exists a matching file name in the chromatin_tracks_input_dir
            motifs_overlapping_tracks_file = params['motifs_overlapping_tracks_output_dir']+'/' + '.'.join(motif_file.split('/')[-1].split('.')[0:-1])+'_overlapping_tracks' + '.bed7'
            scored_motifs_chromatin_tracks_output_file = '.'.join(motifs_overlapping_tracks_file.split('.')[0:-1]) + '_scored.bed10' 
            if not (os.path.exists(motifs_overlapping_tracks_file) and os.path.exists(scored_motifs_chromatin_tracks_output_file)):
                if get_value(params['run_in_parallel_param']) and len(motif_files)>1:
                    p.apply_async(overlay_resources_score_motifs, args=(params['motif_sites_dir']+'/'+motif_file, 
                                                                     params['all_chromatin_makrs_all_cells_combined_dir_path']+'/'+motif_file.split('/')[-1], 
                                                                     scored_motifs_chromatin_tracks_output_file, 
                                                                     motifs_overlapping_tracks_file,
                                                                     normal_expression_per_tissue_origin_per_TF, 
                                                                     matching_cell_name_representative_dict, motifTFName_TFNames_matches_dict, 
                                                                     cells_assays_dict, cell_tfs, tf_cells, assay_cells_datatypes, header))
                else:
                    overlay_resources_score_motifs(params['motif_sites_dir']+'/'+motif_file, 
                                                params['all_chromatin_makrs_all_cells_combined_dir_path']+'/'+motif_file.split('/')[-1], 
                                                scored_motifs_chromatin_tracks_output_file, 
                                                motifs_overlapping_tracks_file,
                                                normal_expression_per_tissue_origin_per_TF,
                                                matching_cell_name_representative_dict, motifTFName_TFNames_matches_dict, 
                                                cells_assays_dict, cell_tfs, tf_cells, assay_cells_datatypes, header)
            motifs_overlapping_tracks_files.append(motifs_overlapping_tracks_file)
            scored_motifs_overlapping_tracks_files.append(scored_motifs_chromatin_tracks_output_file)
    if get_value(params['run_in_parallel_param']) and len(motif_files)>1:
        p.close()
        p.join()
    return motifs_overlapping_tracks_files, scored_motifs_overlapping_tracks_files

def overlay_resources_score_motifs(motif_sites_input_file, chromatin_tracks_input_file, 
                                   scored_motifs_chromatin_tracks_output_file, motifs_overlapping_tracks_file,
                                   normal_expression_per_tissue_origin_per_TF, 
                                   matching_cell_name_representative_dict, motifTFName_TFNames_matches_dict, cells_assays_dict, cell_tfs, tf_cells, assay_cells_datatypes, header): 
    """intersect motifs with chromatin tracks, sort and group the tracks per motif
    Input: moitf instances file (motif pos, name_id, scorePval, strand)
           chromatin data collection file in bed4 format; track pos, track cell#assaytype#value or cell#TFname in case of chip-seq
    Return a file in bed7 format (motif info (6cols), overlapping_tracks. 
    """
    print(("in overlay_resources_score_motifs: " + scored_motifs_chromatin_tracks_output_file))
    
    if not os.path.exists(motifs_overlapping_tracks_file):#intersect motifs and chromatin data
        motifs_chromatin_tracks_output_file_temp_sorted = motifs_overlapping_tracks_file + '_temp_sorted'
        motif_sites_file_obj = BedTool(motif_sites_input_file)
        all_chromatin_makrs_all_cells_file_obj = BedTool(chromatin_tracks_input_file)
        print(("intersecting: " + motif_sites_input_file + ' and ' + chromatin_tracks_input_file))
        motifs_chromatin_tracks_output_file_temp = motif_sites_file_obj.intersect(all_chromatin_makrs_all_cells_file_obj, wo=True)
        os.system('cut -f1-6,10 ' + motifs_chromatin_tracks_output_file_temp.fn + ' | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 -k6,6 > ' + motifs_chromatin_tracks_output_file_temp_sorted)
        motifs_chromatin_tracks_output_file_temp.delete_temporary_history(ask=False)
        motif_sites_input_file_temp_sorted_obj = BedTool(motifs_chromatin_tracks_output_file_temp_sorted)
        k = motif_sites_input_file_temp_sorted_obj.groupby(g=list(range(1,7)), c=7, o=['distinct']).saveas(motifs_overlapping_tracks_file)
        k.delete_temporary_history(ask=False)
        if os.path.exists(motifs_chromatin_tracks_output_file_temp.fn):
            os.remove(motifs_chromatin_tracks_output_file_temp.fn)
        if os.path.exists(motifs_chromatin_tracks_output_file_temp_sorted):
            os.remove(motifs_chromatin_tracks_output_file_temp_sorted)
    
    if not os.path.exists(scored_motifs_chromatin_tracks_output_file):#score each motif-track_overlapping file file
        print(("computing scores to: " + scored_motifs_chromatin_tracks_output_file))
        score_motifs_per_cell(motifs_overlapping_tracks_file, scored_motifs_chromatin_tracks_output_file,
                              normal_expression_per_tissue_origin_per_TF, 
                              matching_cell_name_representative_dict, motifTFName_TFNames_matches_dict, cells_assays_dict, cell_tfs, tf_cells, assay_cells_datatypes, header)
    return motifs_overlapping_tracks_file, scored_motifs_chromatin_tracks_output_file

#Section 3: Score motifs

def score_motifs_per_cell(motifs_overlapping_tracks_file, scored_motifs_chromatin_tracks_output_file,
                          normal_expression_per_tissue_origin_per_TF, 
                          matching_cell_name_representative_dict, motifTFName_TFNames_matches_dict, cells_assays_dict, cell_tfs, tf_cells, assay_cells_datatypes, header,  
                          index_track_names=6, index_motif_name=3): #,run_training = True, weights_per_param_dict = {}, log_base=10, header=True):
    """
    Input: a list of motifs overlapping cell tracks in bed7 format
           normal gene expression dictionary: keys are cell#TF and values are expression levels (float)
           
    Return: list of scored motifs files 
    """
    sep = '\t'
    with open(motifs_overlapping_tracks_file, 'r') as motifs_overlapping_tracks_readfile, open(scored_motifs_chromatin_tracks_output_file, 'w') as scored_motifs_writefile:
        line = motifs_overlapping_tracks_readfile.readline()
        if header:
            header_line = ['posrange', 'chr', 'motifstart', 'motifend', 'name', 'score', 'pval','strand']
            for cell in sorted(cells_assays_dict.keys()):
                for assay in sorted(cells_assays_dict[cell].keys()):
                    header_line.append('_'.join(((cell + "___" + assay).replace('(','').replace(')','')
                                         .replace('-','__')).split()))
            scored_motifs_writefile.write('\t'.join(header_line) + '\n')
        while line:
            split_line = line.strip().split(sep)
            reset_cells_assays_dict = reset_cells_assays_matrix(split_line[index_motif_name].split('_')[0].upper(), cells_assays_dict, cell_tfs, tf_cells, motifTFName_TFNames_matches_dict, assay_cells_datatypes)
            scored_motif_per_cell_per_assay = get_motif_score(split_line, normal_expression_per_tissue_origin_per_TF, 
                          matching_cell_name_representative_dict, motifTFName_TFNames_matches_dict, reset_cells_assays_dict, index_track_names, index_motif_name)  #, run_training, weights_per_param_dict, log_base)
            field_values = process_scored_motif_per_cell_per_assay(split_line[0:index_track_names], scored_motif_per_cell_per_assay)
            scored_motifs_writefile.write('\t'.join(field_values) + '\n')
            line = motifs_overlapping_tracks_readfile.readline()
    return scored_motifs_chromatin_tracks_output_file

def get_motif_score(split_line, normal_expression_per_tissue_origin_per_TF, matching_cell_name_representative_dict, 
                    motifTFName_TFNames_matches_dict, cells_assays_dict, index_track_names, index_motif_name):  #, run_training, weights_per_param_dict, log_base)
    """Calculates a score for a given motif per cell line."""
    #fill in the matrix according to the values in track names column
    tf_name_from_motif = split_line[index_motif_name].split('_')[0].upper()
    
    """Get expression value for the current TF in all tissues"""
    for representative_cell in cells_assays_dict:
        try:
            if 'TFExpr' in cells_assays_dict[representative_cell]:
                if normal_expression_per_tissue_origin_per_TF[representative_cell][tf_name_from_motif]!='NaN':
                    cells_assays_dict[representative_cell]['TFExpr'] = float(normal_expression_per_tissue_origin_per_TF[representative_cell][tf_name_from_motif])
        except ValueError:
            pass
    
    for trackname in split_line[index_track_names].split(','):
        ts = trackname.split('#')
        matching_tissues_cell = []
        #check for matching cell names
        try:
            matching_tissues_cell = matching_cell_name_representative_dict[ts[0]]#HepG2: [Liver,...]
        except KeyError:#skip tracks of cells that have no matching in the rep_cell dict file
            continue
        
        for matching_tissue_cell in matching_tissues_cell:
            if len(ts)==2:
                cells_assays_dict[matching_tissue_cell][ts[1]] = 1
            elif len(ts)==3  and ts[1] != "TFBinding":
                if cells_assays_dict[matching_tissue_cell][ts[1]] == 0.0 or cells_assays_dict[matching_tissue_cell][ts[1]]=='NO':
                    try:
                        cells_assays_dict[matching_tissue_cell][ts[1]] = float(ts[2])
                    except ValueError:
                        cells_assays_dict[matching_tissue_cell][ts[1]] = ts[2]
            elif ts[1] == "TFBinding" and ( len(ts)==3  or  len(ts)==4):
                #a sample motif name is: ZBTB18_MA0698.1 (name_id) only the first is the factor name
                if ts[2].upper()==tf_name_from_motif or ts[2].upper() in motifTFName_TFNames_matches_dict[tf_name_from_motif]:
                    binding_value = 1.0
                    if len(ts)==4:
                        binding_value = float(ts[3])
                    cells_assays_dict[matching_tissue_cell][ts[1]]=binding_value
                else:
                    if cells_assays_dict[matching_tissue_cell]['NumOtherTFBinding'] == 0.0:
                        cells_assays_dict[matching_tissue_cell]['NumOtherTFBinding'] = 1.0
                        cells_assays_dict[matching_tissue_cell]['OtherTFBinding'] = [ts[2]]
                    else:
                        cells_assays_dict[matching_tissue_cell]['NumOtherTFBinding'] += 1.0
                        cells_assays_dict[matching_tissue_cell]['OtherTFBinding'].append(ts[2])
    return cells_assays_dict

def process_scored_motif_per_cell_per_assay(motif_info, scored_motif_per_cell_per_assay):
    "Adds values from the dict to a list and imputate values for NaNs from the other tissues when possible "
    
    field_values = ['[{},{})'.format(motif_info[1], str(int(motif_info[2])+1))]
    field_values.append(motif_info[0].replace("X", '23').replace('Y', '24').replace('M', '25').replace('chr', ''))
    field_values.append(str(int(motif_info[1])))
    field_values.append(str(int(motif_info[2])))
    field_values.append(motif_info[3])
    
    field_values.append(motif_info[4].split('P')[0].strip('S'))
    if 'P' in motif_info[4]:
        field_values.append(motif_info[4].split('P')[1])
    field_values.append(motif_info[5])
    
    '''for inf in motif_info:
        try:
            field_values.append(str(int(inf)))
        except ValueError:
            if "chr" in inf:
                field_values.append(inf.replace("X", '23').replace('Y', '24').replace('M', '25').replace('chr', ''))
            elif inf.startswith('S') and 'P' in inf:
                field_values.append(inf.split('P')[0].strip('S'))
                field_values.append(inf.split('P')[1])
            else:
                field_values.append(inf)
    '''
    processed_cells_assays_dict = {}
    for cell in sorted(cells_assays_dict.keys()):
        processed_cells_assays_dict[cell] = {}
        for assay in sorted(cells_assays_dict[cell].keys()):
            value = ""
            if assay == "OtherTFBinding":
                value = ';'.join(set(scored_motif_per_cell_per_assay[cell][assay]))
            else:
                if scored_motif_per_cell_per_assay[cell][assay]!="NaN":
                    try:
                        value = float(scored_motif_per_cell_per_assay[cell][assay])
                    except ValueError:
                        value = scored_motif_per_cell_per_assay[cell][assay]
                else:
                    value = 'NaN'
            processed_cells_assays_dict[cell][assay] = value
            field_values.append(str(value))
    return field_values

def db_setup(cells_assays_dict, assay_cells_datatypes, motif_cols = ['posrange int4range', 'chr INTEGER', 'motifstart INTEGER', 'motifend INTEGER', 'name text', 'score real', 'pval real', 'strand char(1)'], db_name = 'ActiveMotifs.db', cell_table='motifs', db_user_name='husen', db_host_name='localhost'):
    
    conn = ""
    curs = ""
    try:
        conn = psycopg2.connect("dbname={} user={} host={}".format(db_name, db_user_name, db_host_name))
        curs = conn.cursor()
        curs.execute("SELECT exists(SELECT 1 from pg_catalog.pg_database where datname = %s)", (db_name,))    
        if curs.fetchone()[0]:
            print(("Successfully connected to DB: ", db_name))
        
    except psycopg2.DatabaseError as e:
        print(("Error %s" %e))
        print(("Creating DB: ", db_name))
        con_postgres = psycopg2.connect(dbname='postgres', user=db_user_name, host=db_host_name)
        con_postgres.set_isolation_level(ISOLATION_LEVEL_AUTOCOMMIT)
        curs_postgres = con_postgres.cursor()
        curs_postgres.execute('CREATE DATABASE ' + db_name)
        con_postgres.commit()
        curs_postgres.close()
        con_postgres.close()
        
        conn = psycopg2.connect("dbname={} user={} host={}".format(db_name, db_user_name, db_host_name))
        curs = conn.cursor()
        curs.execute("SELECT exists(SELECT 1 from pg_catalog.pg_database where datname = %s)", (db_name,))    
        if curs.fetchone()[0]:
            print(("Successfully created and connected to DB: ", db_name))
    
    field_names = []
    field_names.extend(motif_cols)
    for cell in sorted(cells_assays_dict.keys()):
        for assay in sorted(cells_assays_dict[cell].keys()):
            data_type = 'text'
            try:
                data_type = assay_cells_datatypes[assay]
            except KeyError:
                pass
            field_names.append('_'.join(((cell + "___" + assay).replace('(','').replace(')','')
                                         .replace('-','__')).split()) + " " + data_type)
    #curs.execute("DROP TABLE IF EXISTS {}".format(cell_table))
    
    create_table_stmt = "CREATE TABLE IF NOT EXISTS {} ({});".format(cell_table, ' ,'.join(field_names))
    curs.execute(create_table_stmt)
    conn.commit()
    conn.close()
    return conn

def open_connection(db_name = 'ActiveMotifs.db', db_user_name='husen', db_host_name='localhost'):
    conn = psycopg2.connect("dbname={} user={} host={}".format(db_name, db_user_name, db_host_name))
    #conn.row_factory = psycopg2.Row
    return conn

def close_connection(conn):
    conn.close()

def get_col_names_from_table(table_name, conn):
    curs = conn.cursor()
    curs.execute("select * FROM {} limit 1".format(table_name))
    return [desc[0] for desc in curs.description]
    
def insert_from_file(i_file, n, db_name, db_user_name, db_host_name, cell_table, header,thread_num=0):
    
    n_processed = 0
    conn = open_connection(db_name, db_user_name, db_host_name)
    curs = conn.cursor()
    with open(i_file, 'r') as ifile:
        print(i_file)
        if header:#read an extra line to skip the first one
            list(islice(ifile, 1))
        while True:
            #next_n_lines = list(islice(ifile,n))
            lines_as_lists = [l.strip().split('\t') for l in list(islice(ifile,n))]
            if not lines_as_lists:
                break
            n_processed+=len(lines_as_lists)#next_n_lines)
            #lines_as_lists = [l.strip().split('\t') for l in next_n_lines]
            try:
                value_marks = ['%s' for i in range(0, len(lines_as_lists[0]))]
                try:
                    curs.executemany('insert into {} values({})'.format(cell_table, ', '.join(value_marks)), lines_as_lists)
                    conn.commit()
                    print(("Thread {} for ({}) has processed: {}".format(thread_num, i_file, n_processed)))
                except:
                    print(("Thread {} coundn't insert to DB ({})".format(thread_num, i_file)))
            except (ValueError, IndexError):
                print(("Empty file" + i_file + "in thread" + str(thread_num)))
    curs.close()
    conn.close()
    return

def insert_into_db(db_name, db_user_name, db_host_name,  
                       cell_table, scored_motifs_overlapping_tracks_files, header):#, dir_to_import, keyword_to_check, header):
    
    if get_value(params['run_in_parallel_param']) and len(scored_motifs_overlapping_tracks_files)>1:
        p = Pool(int(params['number_processes_to_run_in_parallel']))
    thread_num=0
    for i_file in scored_motifs_overlapping_tracks_files:#[f for f in glob.glob('{}/*{}*'.format(dir_to_import, keyword_to_check))]
        if get_value(params['run_in_parallel_param']) and len(scored_motifs_overlapping_tracks_files)>1:
            thread_num+=1
            p.apply_async(insert_from_file, args=[i_file, 1000000, db_name, db_user_name, db_host_name, cell_table, header, thread_num])
        else:
            insert_from_file(i_file, 1000000, db_name, db_user_name, db_host_name, cell_table, header)
    
    if get_value(params['run_in_parallel_param']) and len(scored_motifs_overlapping_tracks_files)>1:
        p.close()
        p.join()
    print(("Data insertion into {} is done".format(cell_table)))
    return

def get_tissue_cell_mappings(cell_assays, assay_names, tissue_cell_mappings_file=params['TissueCellInfo_matches_dict'], key_value_sep='=', values_sep=',', cell_assay_sepe=':',
                             motif_cols = ['posrange', 'chr', 'motifstart', 'motifend', 'name', 'score', 'pval', 'strand']):
    
    tissue_cell_allassays = {}
    tissue_cell_assays = {}
    col_list = []
    
    col_list.extend(motif_cols)
    
    with open(tissue_cell_mappings_file, 'r') as ifile:
        lines = ifile.readlines()
        for line in lines:
            if line.startswith('//') or line.startswith('#') or '=' not in line:
                continue
            sl = line.strip().split(key_value_sep)
            key_value = '_'.join(sl[0].strip().replace('(','').replace(')','').replace('-','__').split())
            if key_value not in list(tissue_cell_assays.keys()):
                tissue_cell_assays[key_value] = {}
                tissue_cell_allassays[key_value] = {}
                for assay_name in assay_names:
                    assay_name = '_'.join(assay_name.strip().replace('(','').replace(')','').replace('-','__').split())
                    tissue_cell_allassays[key_value][assay_name] = 'NaN'
                        
            for s in sl[1].split(values_sep):
                cell = s.split(':')[0]
                if cell in list(cell_assays.keys()):
                    cell_updated_name = '_'.join(cell.replace('(','').replace(')','').replace('-','__').split())
                    if ':' in s:
                        assay = s.split(':')[1]
                        if assay in cell_assays[cell]:
                            assay_updated_name = '_'.join(assay.replace('(','').replace(')','').replace('-','__').split())
                            try:
                                tissue_cell_assays[key_value][assay_updated_name].append((cell_updated_name+'___'+assay_updated_name).encode('ascii','ignore').lower())
                            except KeyError:
                                tissue_cell_assays[key_value][assay_updated_name] = [(cell_updated_name+'___'+assay_updated_name).encode('ascii','ignore').lower()]
                            col_list.append((cell_updated_name+'___'+assay_updated_name).encode('ascii','ignore').lower())
                    else:
                        for assay in cell_assays[cell]:
                            assay_updated_name = '_'.join(assay.replace('(','').replace(')','').replace('-','__').split())
                            try:
                                tissue_cell_assays[key_value][assay_updated_name].append((cell_updated_name+'___'+assay_updated_name).encode('ascii','ignore').lower())
                            except KeyError:
                                tissue_cell_assays[key_value][assay_updated_name] = [(cell_updated_name+'___'+assay_updated_name).encode('ascii','ignore').lower()]
                            col_list.append((cell_updated_name+'___'+assay_updated_name).encode('ascii','ignore').lower())
                                
    return col_list, tissue_cell_assays, tissue_cell_allassays

def db_setup_tissues(tissue_cell_assays, assay_cells_datatypes, motif_cols = ['mid INTEGER'], db_name = 'ActiveMotifs.db', db_user_name='husen', db_host_name='localhost'):
    
    conn = psycopg2.connect("dbname={} user={} host={}".format(db_name, db_user_name, db_host_name))
    curs = conn.cursor()
    tissue_cols = {}
    for tissue in sorted(tissue_cell_assays.keys()):
        fields = []
        fields.extend(motif_cols)
        
        tissue_cols[tissue] = []
        for c in motif_cols:
            tissue_cols[tissue].append(c.split(' ')[0])
        for assay in sorted(tissue_cell_assays[tissue].keys()):
            fields.append(assay + ' ' + assay_cells_datatypes[assay.replace('__', '-')])
            tissue_cols[tissue].append(assay.encode('ascii','ignore'))
        #create a table for each : tissue
        #print 'select count(mid) from {} where mid<= 45000000;'.format(tissue)
        #curs.execute('select count(mid) from {} where mid<= 45000000;'.format(tissue))
        #print curs.fetchone()[0]
        curs.execute("DROP TABLE IF EXISTS {}".format(tissue))
        create_table_stmt = "CREATE TABLE IF NOT EXISTS {} ({});".format(tissue, ' ,'.join(fields))
        print(create_table_stmt)
        curs.execute(create_table_stmt)
    curs.close()
    conn.commit()
    conn.close()
    return tissue_cols

def insert_into_tissues(selected_rows, tissue_cell_assays, tissue_cell_allassays, assay_names,
                           cols_to_write_to,
                           cols_to_write_to_allassays,thread_num,
                           db_name, db_user_name, db_host_name):
    
    print(("Thread {} has started".format(thread_num)))
    conn = open_connection(db_name, db_user_name, db_host_name)
    curs_for_insertion = conn.cursor()
    
    tissues_values = {}
    tissues_fields = {}
    for tissue in sorted(tissue_cell_allassays.keys()):
        tissues_fields[tissue] = ['mid']
        tissues_values[tissue] = []
        for assay in sorted(tissue_cell_allassays[tissue].keys()):
            tissues_fields[tissue].append(assay.encode('ascii','ignore').lower())
        
    #values_to_write = []
    t_process = time.time()
    for row in selected_rows:
        #value_current_row = [row['mid']]
        #reset tissue_cell_allassays for every new row
        for tissue in sorted(tissue_cell_allassays.keys()):
            for assay in sorted(tissue_cell_allassays[tissue].keys()):
                tissue_cell_allassays[tissue][assay]='NaN'
                
        for tissue in sorted(tissue_cell_assays.keys()):
            for assay in sorted(tissue_cell_assays[tissue].keys()):
                values = []
                for col in tissue_cell_assays[tissue][assay]:
                    #col = col.encode('ascii','ignore').lower()
                    #if row[col]!="NaN" and row[col]!=float('nan') and row[col]!='nan':
                    if row[col] not in ["NaN", float('NaN'), float('nan'), 'nan']:
                        try:
                            if not math.isnan(float(row[col])):
                                values.append(float(row[col]))
                        except ValueError:
                            values.append(row[col])              
                
                value = 'NaN'
                if len(values)>0:
                    try:
                        value = sum(values)/float(len(values))
                    except TypeError:
                        value = Counter(values).most_common(1)[0][0]
                    
                tissue_cell_allassays[tissue][assay]=value
                #value_current_row.append(value)
                
        #values_to_write.append(value_current_row)
        #impute missing values
        for assay in assay_names:
            assay = '_'.join(assay.replace('(','').replace(')','').replace('-','__').split())
            tissues_with_NaN_values = []
            tissues_with_values = []
            for tissue in list(tissue_cell_allassays.keys()):
                if tissue_cell_allassays[tissue][assay]=="NaN":
                    tissues_with_NaN_values.append(tissue)
                else:
                    try:
                        tissues_with_values.append(float(tissue_cell_allassays[tissue][assay]))
                    except ValueError:
                        tissues_with_values.append(tissue_cell_allassays[tissue][assay])
            
            if len(tissues_with_NaN_values)>0:
                if len(tissues_with_values)>=4:
                    value = 'NaN'
                    try:
                        value = sum(tissues_with_values)/float(len(tissues_with_values))
                    except TypeError:
                        value = Counter(tissues_with_values).most_common(1)[0][0]
                    if value!='NaN':
                        for tissue in tissues_with_NaN_values:
                            tissue_cell_allassays[tissue][assay] = value
        for tissue in sorted(tissue_cell_allassays.keys()):
            values_selected_row = [row['mid']]
            for assay in sorted(tissue_cell_allassays[tissue].keys()):
                values_selected_row.append(tissue_cell_allassays[tissue][assay])
            tissues_values[tissue].append(values_selected_row)
    print(('t_process (func): ', time.time()-t_process))
    
    #insert all collected values to their perspective tables/tissues
    t_insert = time.time()
    for tissue in sorted(tissues_values.keys()):
        s_chars = ','.join('%s' for i in range(0, len(tissues_fields[tissue])))
        field_names=','.join(tissues_fields[tissue])
        
        #curs_for_insertion.executemany('insert into {table_name} ({field_names}) values ({values})'.format(table_name='temp'+tissue, field_names=field_names, values=s_chars), tissues_values[tissue])
        '''execute_values(cur = curs_for_insertion, 
                       sql = 'insert into {table_name} ({field_names}) values (%s)'.format(table_name=tissue, field_names=','.join(tissues_fields[tissue])), 
                       argslist = tissues_values[tissue],
                       template = None,
                       page_size = 100)
        '''
        t_tissues_values = tuple(tissues_values[tissue])
        dataText = ','.join('('+curs_for_insertion.mogrify(s_chars, row) + ')' for row in t_tissues_values)
        curs_for_insertion.execute('insert into {table_name} ({field_names}) values {values}'.format(table_name=tissue, field_names=field_names, values=dataText)) 
        
    print(('t_insert (func): ', time.time()-t_insert))
    conn.commit()
    del tissues_values
    del tissues_fields
    print(("Thread {} is done".format(thread_num)))
    curs_for_insertion.close()
    conn.close()
    return

def populate_tissue_values(tissue_cell_assays, tissue_cell_allassays, assay_names, col_list, table_from,
                           cols_to_write_to = [], cols_to_write_to_allassays = [],
                           number_of_rows_to_load = 100000,
                           db_name = "regmotifs", db_user_name='husen', db_host_name='localhost'):
    
    for tissue in sorted(tissue_cell_assays.keys()):
        for assay in sorted(tissue_cell_assays[tissue].keys()):
            cols_to_write_to.append(tissue+'___'+assay)
    
    for tissue in sorted(tissue_cell_allassays.keys()):
        for assay in sorted(tissue_cell_allassays[tissue].keys()):
            cols_to_write_to_allassays.append(tissue+'___'+assay)
    
    conn = open_connection(db_name, db_user_name, db_host_name)
    #curs_for_count = conn.cursor(name = "countcurs", cursor_factory=DictCursor)
    thread_num = 0
    #curs_for_count.execute('select count(posrange) from {}'.format(table_from))
    num_rows = 85459976#int(curs_for_count.fetchone()[0])#85459976
    #curs_for_count.close()
    print(('total number of rows to be inserted: ', num_rows))
    curs_for_selection = conn.cursor(name = "selectioncurs", cursor_factory=DictCursor)
    curs_for_selection.itersize = number_of_rows_to_load
    t_for_select =  time.time()
    curs_for_selection.execute('select {} from {}'.format(','.join(col_list), table_from))
    print(('t_to_select: ', time.time() - t_for_select))
    t_for_fetch = time.time()
    selected_rows = curs_for_selection.fetchmany(size=number_of_rows_to_load)
    print(('t_to_fetch: ', time.time()-t_for_fetch))
    print(("Selected {} rows for insertion.".format(str(number_of_rows_to_load))))
    t_jobset = time.time()
    if get_value(params['run_in_parallel_param']) and len(scored_motifs_overlapping_tracks_files)>1:
        print('Running in parallel')
        i = 0
        num_cores = int(params['number_processes_to_run_in_parallel'])
        p = Pool(int(params['number_processes_to_run_in_parallel']))
        while i<num_cores:
            p.apply_async(insert_into_tissues, args=(selected_rows, tissue_cell_assays, tissue_cell_allassays, assay_names,
                                   cols_to_write_to, cols_to_write_to_allassays, thread_num,
                                   db_name, db_user_name, db_host_name))
            num_rows -=len(selected_rows)
            print(('num_rows remaining: ', num_rows))
            
            selected_rows = []
            selected_rows = curs_for_selection.fetchmany(size=number_of_rows_to_load)
            
            if len(selected_rows)<=0:
                p.close()
                p.join()
                break
            if i==num_cores-1:
                p.close()
                p.join()
                print(('t_jobset: ', time.time()-t_jobset))
                t_jobset = time.time()
                p = Pool(int(params['number_processes_to_run_in_parallel']))
                i=0
            i+=1
            thread_num+=1
    else:
        print('Running sequentially')
        while selected_rows:
            t_to_insert = time.time()
            insert_into_tissues(selected_rows, tissue_cell_assays, tissue_cell_allassays, assay_names,
                                   cols_to_write_to, cols_to_write_to_allassays, thread_num,
                                   db_name, db_user_name, db_host_name)
            
            print(('t_to_insert: ', time.time()-t_to_insert))
            num_rows -=len(selected_rows)
            print(('num_rows remaining: ', num_rows))
            t_for_fetch = time.time()
            selected_rows = []
            selected_rows = curs_for_selection.fetchmany(size=number_of_rows_to_load)
            thread_num+=1
            print(('t_to_fetch: ', time.time()-t_for_fetch))
            if num_rows<=0:
                print(("All rows are processed and inserted from {} into tissue tables".format(table_from)))
                sys.exit(0)
    
    if num_rows==0:
        print(("All rows are processed and inserted from {} into tissue tables".format(table_from)))
    else:
        print(("Warning!!! There are {} remaining rows to be inserted".format(num_rows)))
    curs_for_selection.close()
    conn.close()
    return
    
def update_table(conn, cell_table, update_chr_col = True, chr_col='chr', update_posrange_col = True, posrange_col = 'posrange', motifstart_col='motifstart', motifend_col='motifend', pval_score_update=True, score_col_name='score', pval_col_name = 'pval'):
    curs = conn.cursor()
    curs.execute("alter table {} add column IF NOT EXISTS {} int4range".format(cell_table, posrange_col))
    
    if update_chr_col:
        curs.execute("update {0} set {1} = replace(replace(replace(replace({1},'chr', ''), 'X', '23'), 'Y', '24'), 'M', '25');".format(cell_table, chr_col))
        curs.execute("alter table {0} alter column {1} SET DATA TYPE integer USING {1}::integer".format(cell_table, chr_col))
    
    if update_posrange_col:
        curs.execute("update {0} set {1} = int4range({2}, {3}, '[]')".format(cell_table, posrange_col, motifstart_col, motifend_col))
        curs.execute("alter table {0} alter column {1} SET NOT NULL".format(cell_table, posrange_col))
    
    if pval_score_update:
        curs.execute("alter table {0} add column IF NOT EXISTS {1} real".format(cell_table, pval_col_name))
        curs.execute("update {0} set {1} = split_part({2}, 'P', 2)::real".format(cell_table, pval_col_name, score_col_name))
        
        curs.execute("update {0} set {1} = replace(split_part({1}, 'P', 1), 'S', '')".format(cell_table,score_col_name))
        curs.execute("alter table {0} alter column {1} SET DATA TYPE real USING {1}::real".format(cell_table, score_col_name))
    conn.commit()
    
    return

def create_index(conn, cell_table, index_name='indexposrange', index_method = 'gist', index_cols = 'posrange'):
    curs = conn.cursor()
    curs.execute("DROP INDEX IF EXISTS {}".format(index_name))
    creat_index_stmt= "CREATE INDEX IF NOT EXISTS {} ON {} using {} ({})".format(index_name, cell_table, index_method, index_cols)
    print(creat_index_stmt)
    curs.execute(creat_index_stmt)
    conn.commit()
    #conn.close()
    return

def table_contains_data(conn, table_name):
    
    curs = conn.cursor()
    curs.execute('select chr from {} limit 1'.format(table_name))
    if curs.fetchone() is not None:
        print(('{} contains data'.format(table_name)))
        return True
    else:
        return False
    
def create_table_stmt_parallel(db_name, tissue, tissuecols, tissuemotifsimputed):
    
    conn = open_connection(db_name)
    curs = conn.cursor()
    print(("Loading for", tissue))
    print(tissuecols)
    print(("create table if not exists {0} as (select {1} from {2})".format(tissue, tissuecols, tissuemotifsimputed)))
    curs.execute("create table if not exists {0} as (select {1} from {2})".format(tissue, tissuecols, tissuemotifsimputed))
    curs.execute('alter table {0} add column if not exists mid serial unique references motifs(mid);'.format(tissue))
    conn.commit()
    print(("Created", tissue))
    curs.close()
    conn.close()

def populate_table_per_tissue(tissues_table_name, db_name):
    
    conn = open_connection(db_name)
    tissues_col_names = get_col_names_from_table(tissues_table_name, conn)
    conn.close()
    tissues = {}
    for col in tissues_col_names:
        if '___' in col:
            if col.split('___')[0] not in list(tissues.keys()):
                tissues[col.split('___')[0]] = [col + ' as ' + col.split('___')[1]]
            else:
                tissues[col.split('___')[0]].append(col + ' as ' + col.split('___')[1])
    p = Pool()
    for tissue in tissues:
        #p.apply_async(create_table_stmt_parallel, args= (db_name, tissue, ','.join(tissues[tissue]), 'tissuemotifsimputed'))
        create_table_stmt_parallel(db_name, tissue, ','.join(tissues[tissue]), 'tissuemotifsimputed')   
    p.close()
    p.join()
    print("All tables are added successfully")

def split_motifs_parallel(db_name, motifs_table, chr, motif_cols):
    conn = open_connection(db_name)
    curs = conn.cursor()
    new_table_name = "chr"+str(chr)+"motifs"
    print(new_table_name)
    curs.execute('drop table if exists {0};'.format(new_table_name))
    print(('create table if not exists {0} as (select {1} from {2} where chr={3});'.format(new_table_name, ','.join(motif_cols), motifs_table, chr)))
    curs.execute('create table if not exists {0} as (select {1} from {2} where chr={3});'.format(new_table_name, ','.join(motif_cols), motifs_table, chr))
    curs.execute('create index if not exists {0}mid on {1} using btree(mid);'.format("chr"+str(chr), new_table_name))
    curs.execute('create index if not exists {0}posrange on {1} using gist(posrange);'.format("chr"+str(chr), new_table_name))
    conn.commit()
    curs.close()
    conn.close()

def split_motifs_table_by_chr(motifs_table, motif_cols, db_name):
    
    chr_names = list(range(1,26))
    print((db_name, motifs_table))
    p = Pool()
    for chr in chr_names:
        p.apply_async(split_motifs_parallel, args = (db_name, motifs_table, chr, motif_cols))
        #split_motifs_parallel(db_name, motifs_table, chr)
    p.close()
    p.join()
    print('All tables are created')
    return

#to update motif tables
def update_motif_pos_pertable(db_name, chr_table):

    conn = open_connection(db_name=db_name)
    curs = conn.cursor()
    cmd =  'update {} set posrange=int4range(lower(posrange)+1,upper(posrange)+1), motifstart=motifstart+1,motifend=motifend+1;'.format(chr_table)
    print(cmd)
    curs.execute(cmd)
    conn.commit()
    curs.close()
    conn.close()
    
    return

def update_motif_pos(db_name):
    
    p = Pool(8)
    chr_names = list(range(1,26))
    for chr in chr_names:
        chr_table = "chr"+str(chr)+"motifs"
        p.apply_async(update_motif_pos_pertable, args=(db_name, chr_table))
    p.close()
    p.join()
    print("updated all the tables")
    return

if __name__ == '__main__':
    if len(list(params.keys()))==0:
        sys.exit(0)
    data_dir = collect_all_data()
    motifTFName_TFNames_matches_dict = retreive_TFFamilyName_for_motifNames()
    normal_expression_per_tissue_origin_per_TF = get_expression_level_per_originType_per_TF(motifTFName_TFNames_matches_dict)
    tissues_with_gene_expression = list(normal_expression_per_tissue_origin_per_TF.keys())
    
    representative_cell_name_matchings_dict, matching_cell_name_representative_dict = retreive_key_values_from_dict_file(params['cell_names_matchings_dict'])
    assay_cells, cell_assays, cell_tfs, tf_cells, assay_cells_datatypes = get_assay_cell_info(data_dir=data_dir, matching_rep_cell_names_dict=matching_cell_name_representative_dict, tissues_with_gene_expression=tissues_with_gene_expression, generated_dicts_output_file=params['all_chromatin_makrs_all_cells_combined_dir_path']+"_generated_dicts.txt")
    assay_names=list(assay_cells.keys())
    cells_assays_dict = generate_cells_assays_matrix(cell_assays, cell_names=list(representative_cell_name_matchings_dict.keys()), assay_cells_datatypes=assay_cells_datatypes, tissues_with_gene_expression=tissues_with_gene_expression)
    header = True
    motifs_overlapping_tracks_files, scored_motifs_overlapping_tracks_files = run_overlay_resources_score_motifs( 
                    normal_expression_per_tissue_origin_per_TF,
                    matching_cell_name_representative_dict, motifTFName_TFNames_matches_dict, 
                    cells_assays_dict, cell_tfs, tf_cells, assay_cells_datatypes, header)
    
    #write results to the main cellmotifs table
    db_name = 'regmotifs'#'testcellmotifsdb'
    cell_table='cellmotifs'#'testcellmotifs'
    db_user_name = 'husen'
    db_host_name = 'localhost'
    db_setup(cells_assays_dict, assay_cells_datatypes, motif_cols = ['mid serial unique', 'posrange int4range', 'chr INTEGER', 'motifstart INTEGER', 'motifend INTEGER', 'name text', 'score real', 'pval real', 'strand char(1)'], 
             db_name = db_name, cell_table=cell_table, db_user_name=db_user_name, db_host_name='localhost')
    
    conn = open_connection(db_name)
    if not table_contains_data(conn, cell_table): #that is to avoid writing over an already existing table content
        print(("Inserting data into: ", cell_table))
        insert_into_db(db_name = db_name, db_user_name=db_user_name, db_host_name='localhost', 
                       cell_table=cell_table, scored_motifs_overlapping_tracks_files=scored_motifs_overlapping_tracks_files, header=header)# dir_to_import=params['motifs_overlapping_tracks_output_dir'], keyword_to_check="_scored.bed10", header=header)
        print(("Creating index on: ", cell_table))
        create_index(conn, cell_table, index_name='indexposrange', index_method = 'gist', index_cols = 'posrange')
    close_connection(conn)
    
    process_tissues = True
    #write results to the tissues (based on cell motifs) table
    if process_tissues:
        print('Creating tissues tables')
        col_list, tissue_cell_assays, tissue_cell_allassays = get_tissue_cell_mappings(cell_assays, assay_names, 
                                                                               motif_cols = ['mid', 'posrange', 'chr', 'motifstart', 'motifend', 'name', 'score', 'pval', 'strand'])
        tissue_cols = db_setup_tissues(tissue_cell_allassays, assay_cells_datatypes, motif_cols = ['mid INTEGER'], 
                         db_name = db_name, db_user_name=db_user_name, db_host_name='localhost')
        col_list.append('mid')
        print('Inserting data into tissues tables')
        populate_tissue_values(tissue_cell_assays, tissue_cell_allassays, assay_names, col_list, table_from=cell_table, 
                               number_of_rows_to_load = 50000,
                               db_name = db_name, db_user_name=db_user_name, db_host_name='localhost')
        
        print("Creating index on tissues tables")
        conn = open_connection(db_name, db_user_name, db_host_name)
        for tissue_table in sorted(tissue_cols.keys()):
            
            create_index(conn, tissue_table, index_name='index'+tissue_table+'mid', index_method = 'btree', index_cols = 'mid')
        conn.close()
    
    #split motif table per chr
    split_motifs = False
    if split_motifs:
        motif_cols = ['mid', 'posrange', 'chr', 'motifstart', 'motifend', 'name', 'score', 'pval', 'strand']
        split_motifs_table_by_chr(motifs_table=cell_table, motif_cols=motif_cols, db_name=db_name)
    
    update_motif_positions = False
    if update_motif_positions:
        update_motif_pos(db_name)
        
