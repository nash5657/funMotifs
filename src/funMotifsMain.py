"""
Created on Nov 13, 2016

@author: Husen M. Umer

Score motifs: collects cell-type specific data from several public resources and generates a cell-type specific score
for each motif instance in the human genome Input: TF PWMs, human genome, TF chip-seq resources, DNase1 resources,
ChromHMM labels, Gene expression, CAGE peaks, HIC domains, HIC loops, Replication domains Output: A list of motif
instances with a functionality score per cell type Process: the module has three sections
1) collects and processes data from the provided resources,
2) combines data from the collections and overlays them with motifs,
3) computes a score for each motif instance
"""
import os
import sys
from pybedtools import BedTool, set_tempdir, cleanup
import argparse

import Utilities
import DataProcessing
import ProcessTFMotifs
import MotifAnnotation
import GenerateMotifsTables
import pandas as pd


def parse_args():
    '''Parse command line arguments'''
    print('funMotifsMain')
    parser = argparse.ArgumentParser(description='funMotifs Main script')
    parser.add_argument('--param_file', default='', help='')
    parser.add_argument('--temp_dir', default='', help='')

    args, unknown = parser.parse_known_args()
    return args


if __name__ == '__main__':

    # to run this program add param_file=main_parameters.conf as an argument
    args = parse_args()

    # Get parameters from the sys.argv and the argument file
    params = Utilities.get_params(args.param_file)

    # set the temp dir for bedtools operations
    set_tempdir(args.temp_dir)

    """Section 1: Collect resources"""

    # Combine all data tracks into a bed4 files one per chr, also record assay types
    data_dir = DataProcessing.collect_all_data(params['all_chromatin_makrs_all_cells_combined_dir_path'],
                                               params['data_tracks'])

    # Retrieves the TF family name for each TF name
    # TODO: unittest for function
    motifTFName_TFNames_matches_dict = ProcessTFMotifs.retreive_TFFamilyName_for_motifNames(
        params['TF_family_matches_file'])

    # Given a GTEX file retrieve gene expression from each tissue for each TF name
    # TODO: unittest for function
    normal_expression_per_tissue_origin_per_TF = ProcessTFMotifs.get_expression_level_per_originType_per_TF(
        motifTFName_TFNames_matches_dict,
        normal_gene_expression_inputfile=params['normal_gene_expression_inputfile'],
        origin_gene_expression_values_outputfile=params['normal_gene_expression_inputfile'] + "_perTissue_perTF",
        index_tissues_names_row_start=2,
        index_gene_names_col=1,
        index_gene_values_start=2,
        sep='\t')

    # get tissues with gene expression
    tissues_with_gene_expression = normal_expression_per_tissue_origin_per_TF.keys()

    # returns cell names to consider and their different names as dictionary
    # TODO: unittest for function
    representative_cell_name_matchings_dict, matching_cell_name_representative_dict = Utilities.retreive_key_values_from_dict_file(
        params['cell_names_matchings_dict'],
        key_value_sep='=',
        values_sep=',')

    # get assay cell info
    # TODO: unittest for function
    assay_cells, cell_assays, cell_tfs, tf_cells, assay_cells_datatypes = DataProcessing.get_assay_cell_info(
        data_dir=params['all_chromatin_makrs_all_cells_combined_dir_path'],
        sep='\t',
        matching_rep_cell_names_dict=matching_cell_name_representative_dict,
        generated_dicts_output_file=params['all_chromatin_makrs_all_cells_combined_dir_path'] + "_generated_dicts.txt",
        tissues_with_gene_expression=tissues_with_gene_expression)

    # get assay cell names
    assay_names = assay_cells.keys()

    # Generate a default dict based on the information obtained from the tracks in data_dir
    # TODO: unittest for function
    cells_assays_dict = DataProcessing.generate_cells_assays_matrix(cell_assays,
                                                                    cell_names=representative_cell_name_matchings_dict.keys(),
                                                                    assay_cells_datatypes=assay_cells_datatypes,
                                                                    tissues_with_gene_expression=tissues_with_gene_expression)

    """ Section 2: Overlap between the generated resources and motifs """

    # TODO: remove hard-coded variable
    header = True

    # summarize overlapping motifs and their annotations
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

    """ Section 3: Score motifs """

    '''
    Annotation Scores:
    Collect motifs from the training sets
    Use the cell_table above to annotate the motifs in the same cell lines
    Run a log model to generate coeff for each annotation
    '''





    """ Section 4. DB generation """

    # write results to the main cellmotifs table
    if Utilities.get_value(params['create_database']):
        import DBUtilities, GenerateCellTable, GenerateTissueTables

        run_in_parallel_param = Utilities.get_value(params['run_in_parallel_param'])
        number_processes_to_run_in_parallel = int(params['number_processes_to_run_in_parallel'])
        db_name = params['db_name'].lower()
        db_user_name = params['db_user_name']
        db_host_name = params['db_host_name']
        cell_table = 'cell_table'
        tissue_cell_mappings_file = params['TissueCellInfo_matches_dict']

        motif_cols = ['mid serial unique', 'posrange int4range', 'chr INTEGER', 'motifstart INTEGER',
                      'motifend INTEGER', 'name text', 'score real', 'pval real', 'strand char(1)']
        motif_cols_names = ['mid', 'posrange', 'chr', 'motifstart', 'motifend', 'name', 'score', 'pval', 'strand']

        DBUtilities.create_db(db_name, db_user_name, db_host_name)

        process_cell_table = Utilities.get_value(params['generate_cell_table'])

        if process_cell_table:
            GenerateCellTable.generate_cell_table(db_name,
                                                  cell_table,
                                                  db_user_name,
                                                  db_host_name,
                                                  cells_assays_dict,
                                                  assay_cells_datatypes,
                                                  run_in_parallel_param=run_in_parallel_param,
                                                  number_processes_to_run_in_parallel=number_processes_to_run_in_parallel,
                                                  header=header,
                                                  scored_motifs_overlapping_tracks_files=scored_motifs_overlapping_tracks_files,
                                                  motif_cols=motif_cols,
                                                  motif_cols_names=motif_cols_names,
                                                  cell_index_name='indexposrange',
                                                  cell_index_method='gist',
                                                  cell_index_cols='posrange'
                                                  )

        process_tissues = Utilities.get_value(params['generate_tissue_tables'])
        tissues_fscores_table = "all_tissues"
        # write results to the tissues (based on cell motifs) table
        if process_tissues:
            print('Dropping and Creating tissues tables')
            GenerateTissueTables.generate_tissue_tables(db_name,
                                                        cell_table,
                                                        db_user_name,
                                                        db_host_name,
                                                        tissues_fscores_table,
                                                        assay_cells_datatypes,
                                                        cell_assays,
                                                        assay_names,
                                                        tissue_cell_mappings_file,
                                                        run_in_parallel_param,
                                                        number_processes_to_run_in_parallel,
                                                        scored_motifs_overlapping_tracks_files,
                                                        motif_cols_names=motif_cols_names,
                                                        number_of_rows_to_load=250000,
                                                        annotation_weights_inputfile=params[
                                                            'annotation_weights_inputfile'],
                                                        skip_negative_weights=Utilities.get_value(
                                                            params['skip_negative_weights']),
                                                        generate_tissue_from_db=Utilities.get_value(
                                                            params['generate_tissue_from_db'])
                                                        )

        # split motif table per chr
        new_table_name = "motifs"

        # if process_tissues:
        #    motifs_table = tissues_fscores_table

        if not DBUtilities.table_contains_data(db_name, db_user_name, db_host_name,
                                               new_table_name):

            if process_cell_table:
                motifs_table = cell_table

                # comment out motif_cols_names as not declared in function definition 
                GenerateMotifsTables.create_motifs_table(db_name,
                                                         db_user_name,
                                                         db_host_name,
                                                         motifs_table=motifs_table,
                                                         # motif_cols=motif_cols,
                                                         # use only the simplified list to prevent syntax error 
                                                         motif_cols=motif_cols_names,
                                                         new_table_name=new_table_name)
            else:
                GenerateMotifsTables.create_motifs_table_from_file(db_name, db_user_name,
                                                                   db_host_name,
                                                                   scored_motifs_overlapping_tracks_files,
                                                                   motif_cols=motif_cols,
                                                                   motif_cols_names=motif_cols_names,
                                                                   new_table_name=new_table_name,
                                                                   run_in_parallel_param=run_in_parallel_param,
                                                                   number_processes_to_run_in_parallel=number_processes_to_run_in_parallel)

            GenerateMotifsTables.motif_names_table(db_name, db_user_name, db_host_name,
                                                   motifs_table=new_table_name,
                                                   motif_names_table="motif_names"
                                                   )

            split_motifs = Utilities.get_value(params['generate_motif_tables'])
            if split_motifs:
                GenerateMotifsTables.split_motifs_table_by_chr(db_name, db_user_name, db_host_name,
                                                               motifs_table=new_table_name,
                                                               motif_cols=motif_cols_names,
                                                               chr_names=range(1, 26))

            # get PFM for motifs
            PFM_table_name = "motifs_pfm"
            fre_per_allele_per_motif_dict = Utilities.get_freq_per_motif(motif_PFM_input_file=params['motif_PFM_file'])
            GenerateMotifsTables.generate_PFM_table(fre_per_allele_per_motif_dict, PFM_table_name, db_name,
                                                    db_user_name, db_host_name,
                                                    cols=['name text', 'position integer', 'allele char(1)',
                                                          'freq numeric'],
                                                    cols_names=['name', 'position', 'allele', 'freq'])
