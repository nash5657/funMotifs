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
from pybedtools import set_tempdir

import Utilities
import DataProcessing
import ProcessTFMotifs
import MotifAnnotation

if __name__ == '__main__':
    '''to run this program add param_file=main_parameters.conf as an argument'''
    
    '''Get parameters from the sys.argv and the argument file'''
    params = Utilities.get_params(sys.argv)
    if len(params.keys())==0:
        print "Usage: python regMotifsMain.py param_file==../conf/main_parameters.conf"
        sys.exit(0)
    
    '''set the temp dir for bedtools operations'''
    if not os.path.exists(params['temp_dir']):
        os.makedirs(params['temp_dir'])  
    set_tempdir(params['temp_dir'])
    
    #Section 1: Collect resources
    data_dir = DataProcessing.collect_all_data(params['all_chromatin_makrs_all_cells_combined_dir_path'], params['data_tracks'])
    
    
    motifTFName_TFNames_matches_dict = ProcessTFMotifs.retreive_TFFamilyName_for_motifNames(params['TF_family_matches_file'])
    
    normal_expression_per_tissue_origin_per_TF = ProcessTFMotifs.get_expression_level_per_originType_per_TF(motifTFName_TFNames_matches_dict, 
                                                                                            normal_gene_expression_inputfile=params['normal_gene_expression_inputfile'],
                                                                                            origin_gene_expression_values_outputfile = params['normal_gene_expression_inputfile'] + "_perTissue_perTF", 
                                                                                            index_tissues_names_row_start = 2, 
                                                                                            index_gene_names_col = 1, 
                                                                                            index_gene_values_start=2, 
                                                                                            sep='\t')
    
    tissues_with_gene_expression = normal_expression_per_tissue_origin_per_TF.keys()
    
    representative_cell_name_matchings_dict, matching_cell_name_representative_dict = Utilities.retreive_key_values_from_dict_file(params['cell_names_matchings_dict'],
                                                                                                                         key_value_sep='=', 
                                                                                                                         values_sep=',')
    
    
    assay_cells, cell_assays, cell_tfs, tf_cells, assay_cells_datatypes = DataProcessing.get_assay_cell_info(data_dir = params['all_chromatin_makrs_all_cells_combined_dir_path'], 
                                                                                              sep='\t', 
                                                                                              matching_rep_cell_names_dict=matching_cell_name_representative_dict, 
                                                                                              generated_dicts_output_file=params['all_chromatin_makrs_all_cells_combined_dir_path']+"_generated_dicts.txt", 
                                                                                              tissues_with_gene_expression = tissues_with_gene_expression)
    
    assay_names=assay_cells.keys()
    
    cells_assays_dict = DataProcessing.generate_cells_assays_matrix(cell_assays, 
                                                     cell_names=representative_cell_name_matchings_dict.keys(), 
                                                     assay_cells_datatypes=assay_cells_datatypes, 
                                                     tissues_with_gene_expression=tissues_with_gene_expression)
    header = True
    
    #Section2 Overlap between the generated resources and motifs
    motifs_overlapping_tracks_files, scored_motifs_overlapping_tracks_files = MotifAnnotation.run_overlay_resources_score_motifs( 
                    params['motif_sites_dir'],
                    params['all_chromatin_makrs_all_cells_combined_dir_path'], 
                    params['motifs_overlapping_tracks_output_dir'],
                    params['run_in_parallel_param'],
                    params['number_processes_to_run_in_parallel'],
                    normal_expression_per_tissue_origin_per_TF,
                    matching_cell_name_representative_dict, 
                    motifTFName_TFNames_matches_dict, 
                    cells_assays_dict, 
                    cell_tfs, 
                    tf_cells, 
                    assay_cells_datatypes, 
                    header)
    
    
    #Section 3: Score motifs
    '''
    Annotation Scores:
    Collect motifs from the training sets
    Use the cell_table above to annotate the motifs in the same cell lines
    Run a log model to generate coeff for each annotation
    '''
    
    
    
    #Section 4. DB generation, Next step 29 Sep
    #write results to the main cellmotifs table
    if Utilities.get_value(params['create_database']):
        import DBUtilities, GenerateCellTable, GenerateTissueTables
        run_in_parallel_param = Utilities.get_value(params['run_in_parallel_param'])
        number_processes_to_run_in_parallel = Utilities.get_value(params['number_processes_to_run_in_parallel'])
        db_name = params['db_name']
        db_user_name = params['db_user_name']
        db_host_name = params['db_host_name'] 
        cell_table = 'cell_table'
        tissue_cell_mappings_file = params['TissueCellInfo_matches_dict']
        
        motif_cols = ['mid serial unique', 'posrange int4range', 'chr INTEGER', 'motifstart INTEGER', 'motifend INTEGER', 'name text', 'score real', 'pval real', 'strand char(1)']
        motif_cols_names = ['mid', 'posrange', 'chr', 'motifstart', 'motifend', 'name', 'score', 'pval', 'strand']
                    
        DBUtilities.create_db(db_name, db_user_name, db_host_name)

        GenerateCellTable.generate_cell_table(db_name,
                    cell_table,
                    db_user_name,
                    db_host_name,
                    cells_assays_dict,
                    assay_cells_datatypes,
                    run_in_parallel_param,
                    number_processes_to_run_in_parallel,
                    header,
                    scored_motifs_overlapping_tracks_files,
                    motif_cols = motif_cols,
                    motif_cols_names = motif_cols_names,
                    cell_index_name='indexposrange', 
                    cell_index_method = 'gist', 
                    cell_index_cols = 'posrange',
                    )
        
        process_tissues = True
        #write results to the tissues (based on cell motifs) table
        if process_tissues:
            print 'Creating tissues tables'
            GenerateTissueTables.generate_tissue_tables(db_name,
                    cell_table,
                    db_user_name,
                    db_host_name,
                    assay_cells_datatypes,
                    cell_assays,
                    assay_names,
                    tissue_cell_mappings_file,
                    run_in_parallel_param,
                    number_processes_to_run_in_parallel,
                    scored_motifs_overlapping_tracks_files,
                    motif_cols_names=motif_cols_names,
                    number_of_rows_to_load=50000
                )
        
        #split motif table per chr
        split_motifs = False
        if split_motifs:
            DBUtilities.split_motifs_table_by_chr(motifs_table=cell_table, motif_cols=motif_cols_names, db_name=db_name)
        
        