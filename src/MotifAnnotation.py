'''
Created on 28 Sep 2017

@author: husensofteng
'''

import os
from multiprocessing import Pool
from pybedtools import BedTool, set_tempdir, cleanup
import glob
from itertools import starmap 
from itertools import product

def reset_cells_assays_matrix(tf_name_from_motif_name, 
                              cells_assays_dict, 
                              cell_tfs, 
                              tf_cells, 
                              motifTFName_TFNames_matches_dict, 
                              assay_cells_datatypes):
    
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

def get_motif_score(split_line, 
                    normal_expression_per_tissue_origin_per_TF, 
                    matching_cell_name_representative_dict, 
                    motifTFName_TFNames_matches_dict,
                    cells_assays_dict, 
                    index_track_names, 
                    index_motif_name):  #, run_training, weights_per_param_dict, log_base)
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

def process_scored_motif_per_cell_per_assay(motif_info, 
                                            scored_motif_per_cell_per_assay,
                                            cells_assays_dict):
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

def score_motifs_per_cell(motifs_overlapping_tracks_file, 
                          scored_motifs_chromatin_tracks_output_file,
                          normal_expression_per_tissue_origin_per_TF, 
                          matching_cell_name_representative_dict, 
                          motifTFName_TFNames_matches_dict, 
                          cells_assays_dict, 
                          cell_tfs, 
                          tf_cells, 
                          assay_cells_datatypes, 
                          header,  
                          index_track_names, 
                          index_motif_name): #,run_training = True, weights_per_param_dict = {}, log_base=10, header=True):
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
            reset_cells_assays_dict = reset_cells_assays_matrix(split_line[index_motif_name].split('_')[0].upper(), 
                                                                cells_assays_dict, 
                                                                cell_tfs, 
                                                                tf_cells, 
                                                                motifTFName_TFNames_matches_dict, 
                                                                assay_cells_datatypes)
            
            scored_motif_per_cell_per_assay = get_motif_score(split_line, 
                                                              normal_expression_per_tissue_origin_per_TF, 
                                                              matching_cell_name_representative_dict, 
                                                              motifTFName_TFNames_matches_dict, 
                                                              reset_cells_assays_dict, 
                                                              index_track_names, 
                                                              index_motif_name)  #, run_training, weights_per_param_dict, log_base)
            
            field_values = process_scored_motif_per_cell_per_assay(split_line[0:index_track_names], 
                                                                   scored_motif_per_cell_per_assay,
                                                                   cells_assays_dict)
            
            scored_motifs_writefile.write('\t'.join(field_values) + '\n')
            line = motifs_overlapping_tracks_readfile.readline()
    return scored_motifs_chromatin_tracks_output_file

def intersect_motif_and_chromatin_marks (chromatin_tracks_input_file_splitted, motif_sites_input_file):
    motifs_chromatin_tracks_output_file_temp_splitted = chromatin_tracks_input_file_splitted + '_intersected'
    motifs_chromatin_tracks_output_file_temp_splitted_tmp = motifs_chromatin_tracks_output_file_temp_splitted + '_tmp'

    all_chromatin_makrs_all_cells_file_obj = BedTool(chromatin_tracks_input_file_splitted)
    motif_sites_file_obj.intersect(all_chromatin_makrs_all_cells_file_obj, wo=True).saveas(motifs_chromatin_tracks_output_file_temp_splitted_tmp)
    os.system('cut -f1-6,10 ' + motifs_chromatin_tracks_output_file_temp_splitted_tmp + ' | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 -k6,6 > ' + motifs_chromatin_tracks_output_file_temp_splitted)
    os.remove(motifs_chromatin_tracks_output_file_temp_splitted_tmp)
    
    
    cleanup()
    return motifs_chromatin_tracks_output_file_temp_splitted

def overlay_resources_score_motifs(motif_sites_input_file, 
                                   chromatin_tracks_input_file, 
                                   scored_motifs_chromatin_tracks_output_file, 
                                   motifs_overlapping_tracks_file,
                                   normal_expression_per_tissue_origin_per_TF, 
                                   matching_cell_name_representative_dict, 
                                   motifTFName_TFNames_matches_dict, 
                                   cells_assays_dict, 
                                   cell_tfs, 
                                   tf_cells, 
                                   assay_cells_datatypes, 
                                   header, run_in_parallel_param, number_processes_to_run_in_parallel): 
    """intersect motifs with chromatin tracks, sort and group the tracks per motif
    Input: moitf instances file (motif pos, name_id, scorePval, strand)
           chromatin data collection file in bed4 format; track pos, track cell#assaytype#value or cell#TFname in case of chip-seq
    Return a file in bed7 format (motif info (6cols), overlapping_tracks. 
    """
    print("in overlay_resources_score_motifs: " + scored_motifs_chromatin_tracks_output_file)
    
    if not os.path.exists(motifs_overlapping_tracks_file):#intersect motifs and chromatin data
        motifs_chromatin_tracks_output_file_temp = motifs_overlapping_tracks_file + '_temp'
        motifs_chromatin_tracks_output_file_temp_sorted = motifs_chromatin_tracks_output_file_temp+ '_sorted'
        
        print("intersecting: " + motif_sites_input_file + ' and ' + chromatin_tracks_input_file)
        
        
        
        os.system( """split -l 1000000 {} {}""" ).format(chromatin_tracks_input_file,chromatin_tracks_input_file+'_tmp')
        chromatin_tracks_input_file_splitted = glob.glob(chromatin_tracks_input_file+'_tmp*')
        
        
        if run_in_parallel_param and len(motif_files)>1:
            motif_sites_file_obj = BedTool(motif_sites_input_file)
            
            
            
            pm = Pool(number_processes_to_run_in_parallel)
            motifs_chromatin_tracks_output_file_temp_files = pm.starmap(intersect_motif_and_chromatin_marks, product(chromatin_tracks_input_file_splitted, [motif_sites_file_obj]))
            pm.close()
            pm.join()
            
            
   
            
        

    
    

        if not os.path.exists(motifs_chromatin_tracks_output_file_temp):
            with open(motifs_chromatin_tracks_output_file_temp, 'w') as motifs_chromatin_tracks_output_file_temp_outfile:
                for motifs_chromatin_tracks_output_file_temp_file in motifs_chromatin_tracks_output_file_temp_files:
                    with open(motifs_chromatin_tracks_output_file_temp_file, 'r') as motifs_chromatin_tracks_output_file_temp_ifile:
                        motifs_chromatin_tracks_output_file_temp_outfile.write(motifs_chromatin_tracks_output_file_temp_ifile.read())
                       # os.remove(motifs_chromatin_tracks_output_file_temp_file)
        

        os.system( ' sort -k1,1 -k2,2n -k3,3n' + motifs_chromatin_tracks_output_file_temp +'> ' + motifs_chromatin_tracks_output_file_temp_sorted)
        motifs_chromatin_tracks_output_file_temp_fn.delete_temporary_history(ask=False)
        if os.path.exists(motifs_chromatin_tracks_output_file_temp):
            os.remove(motifs_chromatin_tracks_output_file_temp)
        motif_sites_input_file_temp_sorted_obj = BedTool(motifs_chromatin_tracks_output_file_temp_sorted)
        k = motif_sites_input_file_temp_sorted_obj.groupby(g=[1,2,3,4,5,6,7], c=7, o=['distinct']).saveas(motifs_overlapping_tracks_file)
        k.delete_temporary_history(ask=False)
        
        if os.path.exists(motifs_chromatin_tracks_output_file_temp_sorted):
            os.remove(motifs_chromatin_tracks_output_file_temp_sorted)
    
    if not os.path.exists(scored_motifs_chromatin_tracks_output_file):#score each motif-track_overlapping file file
        print("computing scores to: " + scored_motifs_chromatin_tracks_output_file)
        score_motifs_per_cell(motifs_overlapping_tracks_file, 
                              scored_motifs_chromatin_tracks_output_file,
                              normal_expression_per_tissue_origin_per_TF, 
                              matching_cell_name_representative_dict, 
                              motifTFName_TFNames_matches_dict, 
                              cells_assays_dict, 
                              cell_tfs, 
                              tf_cells, 
                              assay_cells_datatypes, 
                              header,
                              index_track_names=6, 
                              index_motif_name=3)
    cleanup()   
    return motifs_overlapping_tracks_file, scored_motifs_chromatin_tracks_output_file


def run_overlay_resources_score_motifs(motif_sites_dir,
                                       all_chromatin_makrs_all_cells_combined_dir_path, 
                                       motifs_overlapping_tracks_output_dir,
                                       run_in_parallel_param,
                                       number_processes_to_run_in_parallel,
                                       normal_expression_per_tissue_origin_per_TF,
                                       matching_cell_name_representative_dict, 
                                       motifTFName_TFNames_matches_dict, 
                                       cells_assays_dict, 
                                       cell_tfs, 
                                       tf_cells, 
                                       assay_cells_datatypes, 
                                       header):
    """pairs matching chromosomes in motif_sites_input_dir and all_chromatin_makrs_all_cells_input_dir and calls overlay_resources_score_motifs
    Input: moitf instances input dir (one file per chr)
           chromatin data collection dir (one file per chr, bed4 format; track pos, track cell#assaytype#value or cell#TFname in case of chip-seq) 
    Return: a list of motif_overlapping_track files
    Precondition: files in motif_sites_input_dir and chromatin_tracks_input_dir should have the same names 
                  Recommended: name files in both dirs as chrNumber, chrX or chrY (where number is between 1-22)
    """
    motif_files = []
    if not os.path.isdir(motif_sites_dir) and os.path.isfile(motif_sites_dir):
        motif_files = [motif_sites_dir]
        motif_sites_dir = "."
    else:
        motif_files = os.listdir(motif_sites_dir)
    
    chromatin_tracks_files = os.listdir(all_chromatin_makrs_all_cells_combined_dir_path)
    if not os.path.exists(motifs_overlapping_tracks_output_dir):
        os.makedirs(motifs_overlapping_tracks_output_dir)
    motifs_overlapping_tracks_files = []
    scored_motifs_overlapping_tracks_files = []
    if run_in_parallel_param and len(motif_files)>1:
        p = Pool(int(number_processes_to_run_in_parallel))
    for motif_file in motif_files:
        
        chr_n_file = motif_file.split('/')[-1]
        with open(motif_sites_dir+'/'+motif_file) as f:
            chr_n_file = f.readline().strip().split('\t')[0].strip()+'.bed'
        if (chr_n_file in chromatin_tracks_files):#it is assumed for every motif file name there exists a matching file name in the chromatin_tracks_input_dir
            motifs_overlapping_tracks_file = motifs_overlapping_tracks_output_dir+'/' + '.'.join(motif_file.split('/')[-1].split('.')[0:-1])+'_overlapping_tracks' + '.bed7'
            scored_motifs_chromatin_tracks_output_file = '.'.join(motifs_overlapping_tracks_file.split('.')[0:-1]) + '_scored.bed10' 
            if not (os.path.exists(motifs_overlapping_tracks_file) and os.path.exists(scored_motifs_chromatin_tracks_output_file)):
                #if run_in_parallel_param and len(motif_files)>1:
                #    p.apply_async(overlay_resources_score_motifs, args=(motif_sites_dir+'/'+motif_file, 
                #                                                     all_chromatin_makrs_all_cells_combined_dir_path+'/'+chr_n_file, 
                #                                                     scored_motifs_chromatin_tracks_output_file, 
                #                                                     motifs_overlapping_tracks_file,
                #                                                     normal_expression_per_tissue_origin_per_TF, 
                #                                                     matching_cell_name_representative_dict, motifTFName_TFNames_matches_dict, 
                #                                                     cells_assays_dict, cell_tfs, tf_cells, assay_cells_datatypes, header))
                #else:
                    overlay_resources_score_motifs(motif_sites_dir+'/'+motif_file, 
                                                all_chromatin_makrs_all_cells_combined_dir_path+'/'+chr_n_file, 
                                                scored_motifs_chromatin_tracks_output_file, 
                                                motifs_overlapping_tracks_file,
                                                normal_expression_per_tissue_origin_per_TF,
                                                matching_cell_name_representative_dict, motifTFName_TFNames_matches_dict, 
                                                cells_assays_dict, cell_tfs, tf_cells, assay_cells_datatypes, header,run_in_parallel_param, number_processes_to_run_in_parallel)
            motifs_overlapping_tracks_files.append(motifs_overlapping_tracks_file)
            scored_motifs_overlapping_tracks_files.append(scored_motifs_chromatin_tracks_output_file)
    if run_in_parallel_param and len(motif_files)>1:
        p.close()
        p.join()
    return motifs_overlapping_tracks_files, scored_motifs_overlapping_tracks_files
    
