"""
Created on Sep 26, 2022

@author: Mark Melzer

Test functions to check section 1 in funMotifsMain.py
"""

import argparse
import unittest
import difflib
import sys

sys.path.insert(1, '../src/')
import DataProcessing
import ProcessTFMotifs
import Utilities


class TestSection1(unittest.TestCase):
    '''def test_collect_all_data1(self):
        """ Test with existing data directory """
        args = parse_args()
        params = Utilities.get_params(args.param_file)
        data_dir = DataProcessing.collect_all_data(params['all_chromatin_makrs_all_cells_combined_dir_path'],
                                                   params['data_tracks'])

        return'''

    '''
    def test_collect_all_data2(self):
        """ Test without existing data directory """
        data_tracks = "./InputTestFilesSection1/DataTracks/CAGE_expr_per_peak_all_cells_promoters_enhancers.bed4,./InputTestFilesSection1/DataTracks/RoaDomainsAllGrouped.bed4,./InputTestFilesSection1/DataTracks/RoaLoopsAllGrouped.bed4,./InputTestFilesSection1/DataTracks/ReplicationDomains.bed4,./InputTestFilesSection1/DataTracks/*ChIP-seq.bed4,./InputTestFilesSection1/DataTracks/*_DNase-seq.bed4,./InputTestFilesSection1/DataTracks/*_ChromatinStates.bed4"
        all_chromatin_makrs_all_cells_combined_dir_path2='./ InputTestFilesSection1 / chromatin_marks_all_cells_onlynarrowpeaks2'
        data_dir = DataProcessing.collect_all_data(all_chromatin_makrs_all_cells_combined_dir_path2, data_tracks)
        # check if the created file is the expected output

        with open('InputTestFilesSection1/chromatin_marks_all_cells_onlynarrowpeaks/chr10.bed', 'r') as a, open(data_dir + '/chr10.bed', 'r') as b:
            differ = difflib.Differ()
            for line in differ.compare(a.readlines(), b.readlines()):
                print(line)
                self.assertNotEqual(line[0], '-')
                self.assertNotEqual(line[0], '+')
        return
    '''

    def test_retreive_TFFamilyName_for_motifNames(self):
        # TODO: create further working and not working tests (space-separated, self written column?)
        outcome = {'KEY1': ['KEY1', 'VALUE1'], 'KEY2': ['KEY2', 'VALUE2A', 'VALUE2B'],
                   'KEY3': ['KEY3', 'VALUE3A', 'VALUE3B'], 'KEY4': ['KEY4', 'VALUE4'],
                   'KEY5': ['KEY5', 'VALUE5A', 'VALUE5B']}
        TF_family_matches_file = "./InputTestFilesSection1/TFNames_motifNames_mapping"
        # for unittest purpose when function below is not needed:
        # TF_family_matches_file = "./InputTestFilesSection1/TFNames_motifNames_mapping"
        x = ProcessTFMotifs.retreive_TFFamilyName_for_motifNames(TF_family_matches_file)
        # assert x == outcome
        return x

    def test_get_expression_level_per_originType_per_TF(self):
        inputfile = \
            "./InputTestFilesSection1/GeneExp/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm_test.gct"
        motifTFName_TFNames_matches_dict = self.test_retreive_TFFamilyName_for_motifNames()
        x = ProcessTFMotifs.get_expression_level_per_originType_per_TF(
            motifTFName_TFNames_matches_dict,
            normal_gene_expression_inputfile=inputfile,
            origin_gene_expression_values_outputfile=inputfile + "_perTissue_perTF",
            index_tissues_names_row_start=2,
            index_gene_names_col=1,
            index_gene_values_start=2,
            sep='\t')
        # TODO: create output that actually is output
        return x

    # TODO: write more test files
    def test_retreive_key_values_from_dict_file(self):
        # use this file when testing this function only and use assert statements
        # infile = "./InputTestFilesSection1/cell_names_to_consider_test.txt"
        infile = "./InputTestFilesSection1/cell_names_to_consider.txt"
        # output_dict = {'SWE': ['Sweden', 'Sverige'], 'MEX': ['Mexico'], 'GER': ['Germany']}
        # output_rev = {'SWE': ['SWE'], 'Sweden': ['SWE'], 'Sverige': ['SWE'], 'MEX': ['MEX'], 'Mexico': ['MEX'], 'GER': ['GER'], 'Germany': ['GER']}
        x = Utilities.retreive_key_values_from_dict_file(
            infile,
            key_value_sep='=',
            values_sep=',')
        # assert output_dict == x[0]
        # assert output_rev == x[1]
        return x

    def test_get_assay_cell_info(self):
        data_dir = "./InputTestFilesSection1/chromatin_marks_all_cells_onlynarrowpeaks/"
        matching_cell_name_representative_dict = self.test_retreive_key_values_from_dict_file()[1]
        tissues_with_gene_expression = self.test_get_expression_level_per_originType_per_TF().keys()
        x = DataProcessing.get_assay_cell_info(
            data_dir=data_dir,
            sep='\t',
            matching_rep_cell_names_dict=matching_cell_name_representative_dict,
            generated_dicts_output_file=data_dir + "_generated_dicts.txt",
            tissues_with_gene_expression=tissues_with_gene_expression)
        # for i in x:
            # print(i)
        # print(x)  assay_cells, cell_assays, cell_tfs, tf_cells, assay_cells_datatypes
        # TODO: implement assert statements
        return x

    def test_generate_cells_assays_matrix(self):
        _, cell_assays, _, _, assay_cells_datatypes = self.test_get_assay_cell_info()
        x = DataProcessing.generate_cells_assays_matrix(cell_assays,
                                                        cell_names=self.test_retreive_key_values_from_dict_file()[
                                                            0].keys(),
                                                        assay_cells_datatypes=assay_cells_datatypes,
                                                        tissues_with_gene_expression=self.test_get_expression_level_per_originType_per_TF().keys())
        print(x)
        # TODO: implement assert statements
        return x


if __name__ == '__main__':
    unittest.main()
