"""
Created on Nov 13, 2016

@author: Husen M. Umer
@contributor: Mark Melzer

Score motifs: collects cell-type specific data from several public resources and generates a cell-type specific score
for each motif instance in the human genome Input: TF PWMs, human genome, TF chip-seq resources, DNase1 resources,
ChromHMM labels, Gene expression, CAGE peaks, HIC domains, HIC loops, Replication domains
Output: A list of motif instances with a functionality score per cell type
Process: the module has six sections
1) collects and processes data from the provided resources,
2) combines data from the collections and overlays them with motifs,
3) write Annotations into database
4) computes a score for each motif instance
5) computes the functional motifs per tissue
6) overlaps the functional motifs with non-coding mutations
"""
import argparse
import os
import psycopg2
import json
import re

import numpy as np
from pybedtools import set_tempdir

import DataProcessing
import FindMotifMutations
import GenerateMotifsTables
import GetFunMotifperTissue as gfmt
import MotifAnnotation
import ProcessTFMotifs
import Utilities
import WeightFeatures
import Entropy
import glob


def parse_args():
    """Parse command line arguments"""
    print('funMotifsMain')
    parser = argparse.ArgumentParser(description='Pipeline for funMotifs. Use the flags -M, -R, -F, and -V to specify,'
                                                 'which parts of the pipeline should be executed. If you use none of'
                                                 'those, nothing will be executed. Also, specify a file containing the '
                                                 'necessary parameters and a temporary directory.')
    parser.add_argument('-m', '--annotateMotifs', action='store_true', help="Set this flag to annotate Motifs and save"
                                                                            " them to the database. Otherwise the"
                                                                            " program assumes that a data base of"
                                                                            " annotated motifs exists already.")
    parser.add_argument('-r', '--regression', action='store_true', help="Set this flag to run the logistic regression."
                                                                        " Otherwise previously computed weights will be"
                                                                        " used from the specified directory.")
    parser.add_argument('-f', '--findFunMotifs', action='store_true', help="Set this flag to identify functional motifs"
                                                                           " based on input data and regression "
                                                                           "output.")
    parser.add_argument('-v', '--findFunMotifVariants', action='store_true', help="Set this flag to overlap the"
                                                                                  " functional motifs with the given"
                                                                                  " genomic variants.")
    parser.add_argument('-p', '--param_file', help='', required=True)
    parser.add_argument('-d', '--temp_dir', default='', help='')
    parser.add_argument('--force_overwrite', default=False, type=bool, help="If true, the specified output path will "
                                                                            "be overwritten if it already exists. Else"
                                                                            " the files in the existing path will be"
                                                                            " used in further functions.")

    args, unknown = parser.parse_known_args()
    return args


# Define the function to save parameters to a file
def save_parameters_to_file(params, file_path):
    # Convert parameters to a dictionary if they are not already in one
    # This step depends on the format of your params
    if not isinstance(params, dict):
        params_dict = params.to_dict()
    else:
        params_dict = params

    # Write the dictionary to a file in JSON format
    with open(file_path, 'w') as file:
        json.dump(params_dict, file, indent=4)
    
    print(f"Parameters successfully saved to {file_path}")

#extract values of annotation_wights file to a list to use in section 5
def extract_variables_from_file(file_path):
    """ Extracts variables and their values from a file and returns them as a dictionary. """
    variable_dict = {}

    # Regular expression pattern to identify variable names and values
    pattern = r'([\w/,]+)=([-+]?\d*\.\d+|\d+)'

    with open(file_path, 'r') as file:
        for line in file:
            # Find all matches in the line
            matches = re.findall(pattern, line)
            for vars, val in matches:
                # Splitting multiple variables if they exist
                vars_split = vars.split(',')
                for var in vars_split:
                    # Adding to the dictionary
                    variable_dict[var] = float(val)

    return variable_dict


if __name__ == '__main__':

    # Establish database connection
    database_name = 'funmotifsdb'
    user_name = 'naser'
    host_address = '127.0.0.1'
    conn = psycopg2.connect(database=database_name, user=user_name, host=host_address)
    cursor = conn.cursor()

    # to run this program add param_file=main_parameters.conf as an argument
    args = parse_args()

    # Get parameters from the sys.argv and the argument file
    params = Utilities.get_params(args.param_file)

    # set the temp dir for bedtools operations
    set_tempdir(args.temp_dir)

    # if force_overwrite, delete results directory to compute section 1 and 2 again
    if Utilities.get_value(args.force_overwrite):
        # TODO: use specified results directory
        os.system("""rm -rf ../results""")

    # get run in parallel
    run_in_parallel_param = Utilities.get_value(params['run_in_parallel_param'])

    """Section 1: Collect resources"""
    # TODO: check which parts of section 1 can be made voluntarily with one of the flags

    # Combine all data tracks into a bed4 files one per chr, also record assay types
    data_dir = DataProcessing.collect_all_data(params['all_chromatin_makrs_all_cells_combined_dir_path'],
                                               params['data_tracks'])

    # Retrieves the TF family name for each TF name
    motifTFName_TFNames_matches_dict = ProcessTFMotifs.retreive_TFFamilyName_for_motifNames(
        params['TF_family_matches_file'])

    # Given a GTEX file retrieve gene expression from each tissue for each TF name
    normal_expression_per_tissue_origin_per_TF = ProcessTFMotifs.get_expression_level_per_originType_per_TF(
        motifTFName_TFNames_matches_dict,
        normal_gene_expression_inputfile=params['normal_gene_expression_inputfile'],
        origin_gene_expression_values_outputfile=params['normal_gene_expression_inputfile'] + "_perTissue_perTF",
        index_tissues_names_row_start=2,
        index_gene_names_col=1,
        index_gene_values_start=2,
        sep='\t')

    # get tissues with gene expression
    tissues_with_gene_expression = list(normal_expression_per_tissue_origin_per_TF.keys())

    # returns cell names to consider and their different names as dictionary
    representative_cell_name, matching_cell_name_representative_dict = Utilities.retreive_key_values_from_dict_file(
        params['cell_names_matchings_dict'],
        key_value_sep='=',
        values_sep=',')

    # get assay cell info
    assay_cells, cell_assays, cell_tfs, tf_cells, assay_cells_datatypes = DataProcessing.get_assay_cell_info(
        data_dir=params['all_chromatin_makrs_all_cells_combined_dir_path'],
        sep='\t',
        matching_rep_cell_names_dict=[representative_cell_name, matching_cell_name_representative_dict],
        generated_dicts_output_file=params['all_chromatin_makrs_all_cells_combined_dir_path'] + "_generated_dicts.txt",
        tissues_with_gene_expression=tissues_with_gene_expression)

    # get assay cell names
    assay_names = list(assay_cells.keys())

    # Generate a default dict based on the information obtained from the tracks in data_dir
    # TODO: cell names below only considers the representative names not the alternative ones --> check
    # work if it is not supposed to contain values , but only 0.0 or 'NO'
    cells_assays_dict = DataProcessing.generate_cells_assays_matrix(cell_assays,
                                                                    cell_names=representative_cell_name,
                                                                    assay_cells_datatypes=assay_cells_datatypes,
                                                                    tissues_with_gene_expression=
                                                                    tissues_with_gene_expression)

    # Generate mapping from cell types to tissue
    matching_tissue_to_cell, cell_name_for_tissue = Utilities.cell_to_tissue_matches(
        params['TissueCellInfo_matches_dict'])

    print("Finished Section 1")

    """ Section 2: Overlap between the generated resources and motifs """

    if args.annotateMotifs:

        # summarize overlapping motifs and their annotations
        motifs_overlapping_tracks_files, scored_motifs_overlapping_tracks_files = \
            MotifAnnotation.run_overlay_resources_score_motifs(params['motif_sites_dir'],
                                                               params[
                                                                   'all_chromatin_makrs_all_cells_combined_dir_path'],
                                                               params['motifs_overlapping_tracks_output_dir'],
                                                               run_in_parallel_param,
                                                               params['number_processes_to_run_in_parallel'],
                                                               normal_expression_per_tissue_origin_per_TF,
                                                               # matching_tissue_to_cell,
                                                               matching_cell_name_representative_dict,
                                                               motifTFName_TFNames_matches_dict,
                                                               cells_assays_dict, cell_tfs, tf_cells,
                                                               assay_cells_datatypes)

        """ Section 3. DB generation """

        # write results to the main cell motifs table
        # TODO: that if statement and parameter should be removed as it is already stated with the flag args.m
        if Utilities.get_value(params['create_database']):
            import DBUtilities, GenerateCellTable, GenerateTissueTables

            # load parameters from parameter file
            number_processes_to_run_in_parallel = int(params['number_processes_to_run_in_parallel'])
            db_name = params['db_name'].lower()
            db_user_name = params['db_user_name']
            db_host_name = params['db_host_name']
            # TODO: remove hardcoded variable
            cell_table = 'cell_table'
            tissue_cell_mappings_file = params['TissueCellInfo_matches_dict']
            db_dir = params['db_dir']
            logfile = params['logfile']

            # set data base table column names for motif information and their names
            motif_cols = ['mid serial unique', 'posrange int4range', 'chr INTEGER', 'motifstart INTEGER',
                          'motifend INTEGER', 'name text', 'score real', 'pval real', 'strand char(1)']
            motif_cols_names = ['mid', 'posrange', 'chr', 'motifstart', 'motifend', 'name', 'score', 'pval', 'strand']

            # TODO: start psql server
            # start the psql server
            # DBUtilities.start_psql_server(db_dir, logfile)
            DBUtilities.create_db(db_name, db_user_name, db_host_name)

            # generate cell tables
            if Utilities.get_value(params['generate_cell_tables']):
                GenerateCellTable.generate_cell_table(db_name,
                                                      cell_table,
                                                      db_user_name,
                                                      db_host_name,
                                                      cells_assays_dict,
                                                      assay_cells_datatypes,
                                                      run_in_parallel_param=run_in_parallel_param,
                                                      number_processes_to_run_in_parallel=number_processes_to_run_in_parallel,
                                                      scored_motifs_overlapping_tracks_files=scored_motifs_overlapping_tracks_files,
                                                      motif_cols=motif_cols,
                                                      motif_cols_names=motif_cols_names,
                                                      cell_index_name='indexposrange',
                                                      cell_index_method='gist',
                                                      cell_index_cols='posrange'
                                                      )

            # TODO: remove hard-coded variable
            tissues_fscores_table = "all_tissues"
            # generate tissue tables
            if Utilities.get_value(params['generate_tissue_tables']):
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

            # TODO: remove hard-coded variable
            # split motif table per chr
            new_table_name = "motifs"

            if not DBUtilities.table_contains_data(db_name, db_user_name, db_host_name,
                                                   new_table_name):

                if Utilities.get_value(params['generate_cell_tables']):
                    motifs_table = cell_table

                    GenerateMotifsTables.create_motifs_table(db_name,
                                                             db_user_name,
                                                             db_host_name,
                                                             motifs_table=motifs_table,
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

                if Utilities.get_value(params['generate_motif_tables']):
                    GenerateMotifsTables.split_motifs_table_by_chr(db_name, db_user_name, db_host_name,
                                                                   motifs_table=new_table_name,
                                                                   motif_cols=motif_cols_names,
                                                                   chr_names=list(range(1, 26)))

                # get PFM for motifs
                PFM_table_name = "motifs_pfm"
                fre_per_allele_per_motif_dict = Utilities.get_freq_per_motif(
                    motif_PFM_input_file=params['motif_PFM_file'])
                GenerateMotifsTables.generate_PFM_table(fre_per_allele_per_motif_dict, PFM_table_name, db_name,
                                                        db_user_name, db_host_name,
                                                        cols=['name text', 'position integer', 'allele char(1)',
                                                              'freq numeric'],
                                                        cols_names=['name', 'position', 'allele', 'freq'])
    else:
        print("Use existing motif annotations from specified data base")
        db_name = params['db_name'].lower()
        db_user_name = params['db_user_name']
        db_host_name = params['db_host_name']
        # TODO: how to integrate database in here, and something else to do in else statement?

    # TODO: check at which point created files can be deleted

    """ Section 4: Score motifs """

    if args.regression:
        cell_table = 'cell_table'
        datafiles_motifs_dir = params['motif_sites_dir']

        training_data_dir = params['trainings_data_dir']

        motif_info_col_names = ['chr', 'motifstart', 'motifend', 'name', 'score', 'pval', 'strand']

        col_names_to_weight_param = ['ChromHMM'.lower(), 'DNase__seq'.lower(), 'FANTOM'.lower(),
                                     'NumOtherTFBinding'.lower(), 'RepliDomain'.lower(), 'TFBinding'.lower(),
                                     'TFExpr'.lower(), 'score'.lower(), 'footprints'.lower(), 'cCRE'.lower(),
                                     'IndexDHS'.lower(), 'RegElem'.lower()]
        logit_params = WeightFeatures.get_param_weights(
        training_data_dir, col_names_to_weight_param, db_name,
        motif_info_col_names, cell_table, db_user_name,
        cell_name_for_tissue, matching_tissue_to_cell,
        motif_split_chr=datafiles_motifs_dir
    )

    #print(logit_params.summary())
    try:
        exp_params = np.exp(logit_params.params)
        print(exp_params)

        # Call function to save parameters
        save_parameters_to_file(exp_params, 'logit_parameters.json')

    except Exception as e:
        print("Error in computing parameters:", e)
        # Handle other exceptions or errors here


    """ Section 5: Get the functional motifs per tissue """
    if args.findFunMotifs:
    #     # Existing code for getting functional motifs
    #     funMotifs_per_Tissue = gfmt.get_functional_motifs_per_tissue(
    #         params=logit_params.params,
    #         tissues=cell_name_for_tissue.keys(),
    #         cursor=cursor  # Passing the cursor instead of db credentials
    #     )
    #     # TODO: save output to database
    # else:
    #     print("Use existing functional Motifs")
    #     # Load functional motifs from file
        try:
            file_path = 'annotation_wights.txt'
            variables = extract_variables_from_file(file_path)
            loaded_params = list(variables.values())
            # Use loaded_params as needed
            # Call function using loaded parameters
            funMotifs_per_Tissue = gfmt.get_functional_motifs_per_tissue(
                params=loaded_params,
                tissues=cell_name_for_tissue.keys(),
                cursor=cursor
        )
        except FileNotFoundError:
            print("Parameters file not found. Please check the file path.")
        except json.JSONDecodeError:
            print("Error decoding the JSON file. Please check the file contents.")

    # Close the cursor and connection after all operations
    cursor.close()
    conn.close()
    

    """ Section 6: Overlap motifs with non-coding mutations and compute entropy """

    if args.findFunMotifVariants:
        # crate output files containing functional motifs containing variant for each tissue
        FindMotifMutations.find_funMotif_variants(funMotifs=funMotifs_per_Tissue, tissues=cell_name_for_tissue.keys(),
                                                  variant_file="",
                                                  output_file="../results/funMotif_variants/overlap_in_",
                                                  db_user_name=db_user_name)
        # compute entropy for each tissue
        for file in glob.glob("../results/funMotif_variants/overlap_in_*"):
            Entropy.compute_entropy(infile=file, outfile=file + "_with_Entropy", db_name=db_name,
                                    db_user_name=db_user_name)
        # TODO: save output to database
    # no else needed, only previous functions were interesting to execute
