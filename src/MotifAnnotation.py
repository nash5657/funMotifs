'''
Created on 28 Sep 2017

@author: husensofteng
'''

import sys, os
from multiprocessing import Pool
from pybedtools import BedTool, set_tempdir, cleanup
import glob
from itertools import starmap 
from itertools import product
from collections import Counter

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
                if tf_name_from_motif in normal_expression_per_tissue_origin_per_TF[representative_cell]:
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
    
    if 'P' in motif_info[4]: 
        field_values.append(motif_info[4].split('P')[0].strip('S'))
        if 'P' in motif_info[4]:
            field_values.append(motif_info[4].split('P')[1])
        field_values.append(motif_info[5])
    else:
        #score
        field_values.append(motif_info[4])
        #p-value
        field_values.append(motif_info[5])
        #strand
        field_values.append(motif_info[6])
        
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
                          normal_expression_per_tissue_origin_per_TF, 
                          matching_cell_name_representative_dict, 
                          motifTFName_TFNames_matches_dict, 
                          cells_assays_dict, 
                          cell_tfs, 
                          tf_cells, 
                          assay_cells_datatypes, 
                          index_track_names, 
                          index_motif_name): #,run_training = True, weights_per_param_dict = {}, log_base=10, header=True):
    """
    Input: a list of motifs overlapping cell tracks in bed7 format
           normal gene expression dictionary: keys are cell#TF and values are expression levels (float)
           
    Return: list of scored motifs files 
    """
    scored_motifs_chromatin_tracks_output_file = motifs_overlapping_tracks_file + '_scored'
    if not os.path.exists(scored_motifs_chromatin_tracks_output_file):
        sep = '\t'
        with open(motifs_overlapping_tracks_file, 'r') as motifs_overlapping_tracks_readfile, open(scored_motifs_chromatin_tracks_output_file, 'w') as scored_motifs_writefile:
            line = motifs_overlapping_tracks_readfile.readline()
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


def overlay_resources_score_motifs(motif_sites_input_file, 
                                   motifs_overlapping_tracks_output_dir,
                                   chromatin_tracks_dir_path,  
                                   chromatin_tracks_files): 
    
    

    """intersect motifs with chromatin tracks, sort and group the tracks per motif
    Input: moitf instances file (motif pos, name_id, scorePval, strand)
           chromatin data collection file in bed4 format; track pos, track cell#assaytype#value or cell#TFname in case of chip-seq
    Return a file in bed7 format (motif info (6cols), overlapping_tracks. 
    """
    
    #for motif_sites_input_file in motif_sites_input_files:
    print("Called here, motifanno, line 231")
    with open(motif_sites_input_file) as f:
        chr_n_file = f.readline().strip().split('\t')[0].strip()+'.bed'
        if (chr_n_file in chromatin_tracks_files):#it is assumed for every motif file name there exists a matching file name in the chromatin_tracks_input_dir
            motifs_overlapping_tracks_file = motifs_overlapping_tracks_output_dir+'/' + '.'.join(motif_sites_input_file.split('/')[-1].split('.')[0:-1])+'_overlapping_tracks' + '.bed7'
            motifs_overlapping_tracks_file_tmp = motifs_overlapping_tracks_file + '_tmp'
            print("in overlay_resources_score_motifs: " + motifs_overlapping_tracks_file)
            print(os.path.exists(motifs_overlapping_tracks_file))
            if not os.path.exists(motifs_overlapping_tracks_file):
                motif_sites_input_file_sorted = motif_sites_input_file + '_sorted'
                chromatin_tracks_input_file = chromatin_tracks_dir_path +'/'+ chr_n_file
                chromatin_tracks_input_file_sorted = chromatin_tracks_input_file + '_sorted'
                
                print("intersecting: " + motif_sites_input_file + ' and ' + chromatin_tracks_input_file)
                
                os.system("""sort -k1,1 -k2,2n -k3,3n {} > {}""".format(motif_sites_input_file, motif_sites_input_file_sorted))
                os.system("""sort -k1,1 -k2,2n -k3,3n {} > {}""".format(chromatin_tracks_input_file, chromatin_tracks_input_file_sorted))
                

                motif_sites_file_obj = BedTool(motif_sites_input_file_sorted)
                motif_sites_file_obj.map(BedTool(chromatin_tracks_input_file_sorted), c=4, o=['collapse']).saveas(motifs_overlapping_tracks_file_tmp)
                
                with open(motifs_overlapping_tracks_file_tmp, 'r') as infile, open(motifs_overlapping_tracks_file, 'w') as outfile:
                        line = infile.readline()
                        # print("Here: ", line)
                        while line:
                            
                            sline = line.split('\t')
                            print("Here, Motifanno, line 258: ", len(sline))
                            if(len(sline)>6):
                                print("here, motifanno, line 260", len(sline))
                                print(sline[7], '\t', sline[8])
                                if(sline[7]!='.'):
                                    my_list=sline[7].split(',')
                                    cell_assay_values_dict_ChromHMM = {}
                                    cell_assay_values_dict_cCRE = {}
                                    cell_assay_values_dict_IndexDHS = {}
                                    cell_assay_values_dict_RegElem = {}
                                    cell_assay_values_dict_DNaseq = {}
                                    elem_list =[]
                                    #elem_list_EpiMap =[]
                                    for elem in my_list:
                                        #print(elem)
                    
                                        cell_value=elem.split('#')[0]
                                        assay_value = elem.split('#')[1]
                                        if(len(elem.split('#'))>2):
                                            state_value = elem.split('#')[2].rstrip("\n")
                    
                                        if assay_value== "ChromHMM":
                                            if cell_value not in cell_assay_values_dict_ChromHMM.keys():
                                                cell_assay_values_dict_ChromHMM[cell_value] = []
    
                                            cell_assay_values_dict_ChromHMM[cell_value].append(state_value)
                                        elif assay_value== "cCRE": 
                                            if cell_value not in cell_assay_values_dict_cCRE.keys():
                                                cell_assay_values_dict_cCRE[cell_value] = []
                                            cell_assay_values_dict_cCRE[cell_value].append(state_value)
                                        elif assay_value== "IndexDHS":
                                            if cell_value not in cell_assay_values_dict_IndexDHS.keys():
                                                cell_assay_values_dict_IndexDHS[cell_value] = []
                                            cell_assay_values_dict_IndexDHS[cell_value].append(state_value)
                                        elif assay_value== "RegElem":
                                            if cell_value not in cell_assay_values_dict_RegElem.keys():
                                                cell_assay_values_dict_RegElem[cell_value] = []
                                            cell_assay_values_dict_RegElem[cell_value].append(state_value)
                                        elif assay_value== "DNase-seq":
                                            if cell_value not in cell_assay_values_dict_DNaseq.keys():
                                                cell_assay_values_dict_DNaseq[cell_value] = []
                                            cell_assay_values_dict_DNaseq[cell_value].append(float(state_value))
                                        else:
    
                                            elem_list.append(elem.rstrip("\n"))                               
                                    for cell in cell_assay_values_dict_ChromHMM:
                                        #print(cell)
                                        #print(cell+"#ChromHMM#"+Counter(cell_assay_values_dict_ChromHMM[cell]).most_common(1)[0][0])
                                        elem_list.append(cell+"#ChromHMM#"+Counter(cell_assay_values_dict_ChromHMM[cell]).most_common(1)[0][0])
                                    for cell in cell_assay_values_dict_cCRE.keys():
                                            #print(cell+"#cCRE#"+Counter(cell_assay_values_dict_cCRE[cell_value]).most_common(1)[0][0])
                                        elem_list.append(cell+"#cCRE#"+Counter(cell_assay_values_dict_cCRE[cell]).most_common(1)[0][0])
    
                                    for cell in cell_assay_values_dict_IndexDHS.keys():
                                            #print(cell_assay_values_dict_IndexDHS[cell])
                                        elem_list.append(cell+"#IndexDHS#"+Counter(cell_assay_values_dict_IndexDHS[cell]).most_common(1)[0][0])
                                    for cell in cell_assay_values_dict_RegElem.keys():
                                            #print(cell_assay_values_dict_IndexDHS[cell])
                                        elem_list.append(cell+"#RegElem#"+Counter(cell_assay_values_dict_RegElem[cell]).most_common(1)[0][0])
                                    for cell in cell_assay_values_dict_DNaseq.keys():
                                            #print(cell_assay_values_dict_IndexDHS[cell])
                                        elem_list.append(cell+"#DNase-seq#"+str(max(cell_assay_values_dict_DNaseq[cell])))
                    
                                    outfile.write('\t'.join(sline[0:7])+'\t'+','.join(elem_list)+'\n')
                    
                            line = infile.readline()
                os.remove(motif_sites_input_file_sorted)
                os.remove(chromatin_tracks_input_file_sorted)
                os.remove(motifs_overlapping_tracks_file_tmp)
            
            
        cleanup()   
    return motifs_overlapping_tracks_file


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

   
    motif_files_full_path =  [motif_sites_dir+'/' + s  for s in motif_files]
    print(motif_files_full_path )
    
    chromatin_tracks_files = os.listdir(all_chromatin_makrs_all_cells_combined_dir_path)
    if not os.path.exists(motifs_overlapping_tracks_output_dir):
        os.makedirs(motifs_overlapping_tracks_output_dir)
    #scored_motifs_chromatin_tracks_output_file = '.'.join(motifs_overlapping_tracks_file.split('.')[0:-1]) + '_scored.bed10' 
    #if not (os.path.exists(motifs_overlapping_tracks_file) and os.path.exists(scored_motifs_chromatin_tracks_output_file)):
    print(run_in_parallel_param, motif_files)
    if run_in_parallel_param and len(motif_files)>1:
        p = Pool(int(number_processes_to_run_in_parallel))
        print("Here, MotifAnnot, Line373")
        print(motif_files_full_path, [motifs_overlapping_tracks_output_dir], [all_chromatin_makrs_all_cells_combined_dir_path], [chromatin_tracks_files])
        motifs_overlapping_tracks_files = p.starmap(overlay_resources_score_motifs, product(motif_files_full_path, 
                                                                    [motifs_overlapping_tracks_output_dir],
                                                                     [all_chromatin_makrs_all_cells_combined_dir_path],
                                                                    [chromatin_tracks_files]))
        print("Here, MotifAnnot, Line378")
        p.close()
        p.join()
    else:
        print("Here, Motifanno, Line383")
        print(motif_files_full_path)
        motifs_overlapping_tracks_files = overlay_resources_score_motifs(motif_files_full_path[0], 
                                                motifs_overlapping_tracks_output_dir,
                                                all_chromatin_makrs_all_cells_combined_dir_path, 
                                                chromatin_tracks_files)

    #score intersected track files
    print("Test: ", motifs_overlapping_tracks_files)
    scored_motifs_overlapping_tracks_files =[]
    for motifs_overlapping_tracks_file in motifs_overlapping_tracks_files:
        scored_motifs_chromatin_tracks_output_file = '.'.join(motifs_overlapping_tracks_file.split('.')[0:-1]) + '_scored.bed10' 
        print("Test: ", motifs_overlapping_tracks_file)
        #if statement inserted for debugging
        print("\'")
        if motifs_overlapping_tracks_file[0] is not "\'":
            motifs_overlapping_tracks_file = "\'" + motifs_overlapping_tracks_file + "\'"
        print("Test: ", motifs_overlapping_tracks_file)
        with open(motifs_overlapping_tracks_file) as f:
                    count = sum(1 for _ in f)
        if not os.path.exists(scored_motifs_chromatin_tracks_output_file):#score each motif-track_overlapping file file
            print("computing scores to: " + scored_motifs_chromatin_tracks_output_file)
            index_track_names=7
            index_motif_name=3
            with open(scored_motifs_chromatin_tracks_output_file, 'w') as scored_motifs_writefile:
                if header:
                    header_line = ['posrange', 'chr', 'motifstart', 'motifend', 'name', 'score', 'pval','strand']
                    for cell in sorted(cells_assays_dict.keys()):
                        for assay in sorted(cells_assays_dict[cell].keys()):
                            if cell[0].isdigit():
                            #if cell=='22Rv1' or cell=='8988T':
                                cell='a'+cell
                                
                            #if cell=="Ammon's horn":
                            #    cell="Ammons horn"
                            #if cell=="Peyer's patch":
                            #    cell="Peyers patch"
                            cell_name ='_'.join(((cell + "___" + assay).replace('(','').replace(')','')
                                                 .replace('-','__').replace('.','').replace("'","")).split())
                            header_line.append('"'+cell_name+'"')
                    #print(header_line)
                    scored_motifs_writefile.write('\t'.join(header_line) + '\n')
            if (run_in_parallel_param):
                os.system( """split -l 200000 {} {}""" .format(motifs_overlapping_tracks_file,motifs_overlapping_tracks_file+'_tmp'))
                motifs_overlapping_tracks_file_splitted = glob.glob(motifs_overlapping_tracks_file+'_tmp*')
                p = Pool(int(number_processes_to_run_in_parallel))
                p.starmap(score_motifs_per_cell, product(motifs_overlapping_tracks_file_splitted, 
                                      [normal_expression_per_tissue_origin_per_TF], 
                                      [matching_cell_name_representative_dict], 
                                      [motifTFName_TFNames_matches_dict], 
                                      [cells_assays_dict], 
                                      [cell_tfs], 
                                      [tf_cells], 
                                      [assay_cells_datatypes], 
                                      [index_track_names], 
                                      [index_motif_name]))
                p.close()
                p.join() 
                #remove tmp splitted files
                with open(scored_motifs_chromatin_tracks_output_file, 'a') as scored_motifs_writefile:
                    for f in motifs_overlapping_tracks_file_splitted:
                        print(f+'_scored')
                        with open(f+'_scored', 'r') as f_score_ifile:
                            l = f_score_ifile.readline()
                            while l:
                                scored_motifs_writefile.write(l)
                                l = f_score_ifile.readline()
                            
                            
                        
                        f_score_ifile.close()
                        os.remove(f)
                        os.remove(f+'_scored')   
                scored_motifs_writefile.close()
            else:
                score_motifs_per_cell(motifs_overlapping_tracks_file, 
                                      normal_expression_per_tissue_origin_per_TF, 
                                      matching_cell_name_representative_dict, 
                                      motifTFName_TFNames_matches_dict, 
                                      cells_assays_dict, 
                                      cell_tfs, 
                                      tf_cells, 
                                      assay_cells_datatypes, 
                                      index_track_names, 
                                      index_motif_name)
                
                    
        scored_motifs_overlapping_tracks_files.append(scored_motifs_chromatin_tracks_output_file)   
    print(scored_motifs_overlapping_tracks_files)
    return motifs_overlapping_tracks_files, scored_motifs_overlapping_tracks_files
    
