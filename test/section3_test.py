"""
Created on Sep 26, 2022

@author: Mark Melzer

Test functions to check section 1 in funMotifsMain.py
"""

import unittest
import sys
import pandas as pd
import statsmodels.api as sm
import difflib


sys.path.insert(1, '../src/')

import WeightFeatures
from GetDBData_TrainingSets import get_cell_info_for_motifs, get_cell_info_for_regions


class TestSection3(unittest.TestCase):
    '''
    def test_funMotifs_logit(self):
        mtcars = sm.datasets.get_rdataset("mtcars", "datasets", cache=True).data
        df = pd.DataFrame(mtcars)
        outcome_col = 'vs'
        col_weight = ['disp', 'wt']
        x = WeightFeatures.funMotifs_logit(df[outcome_col], df[col_weight])

        # get values from R results (fastLR() from RcppNumerical package)
        wt_coeff = 2.48242
        disp_coeff = -0.03977
        log_likelihood = -10.92891

        # get values from R result sm.Logit(method=bfgs)
        py_log_likelihood = x.llf
        py_wt_coeff = x.params['wt']
        py_disp_coeff = x.params['disp']

        # check that results are essentially the same
        assert abs(wt_coeff - py_wt_coeff) < 0.001
        assert abs(disp_coeff - py_disp_coeff) < 0.001
        assert abs(log_likelihood - py_log_likelihood) < 0.001

        return x

    def test_get_coeff(self):
        # test if the code in get_coeff works (funMotifs_Logit checked before)
        df = pd.read_csv("InputTestFilesSection3/combineddfs.tsv_raw.tsv", sep='\t')
        outcome_col = 'activity_score'
        col_names_to_weight_param = df.columns
        col_names_to_weight_param = list(col_names_to_weight_param[3:])
        print(col_names_to_weight_param)


        dfout_filename = "./InputTestFilesSection3/get_coeff_outfile.tsv"

        print('here')
        x = WeightFeatures.get_coeff(df=df, cols_to_weight=df[col_names_to_weight_param], outcome_col=outcome_col,
                                     col_names_to_weight_param=col_names_to_weight_param,
                                     dfout_filename=dfout_filename)
        print('here2')

        assert self.compare_files(dfout_filename, "InputTestFilesSection3/combineddfs.tsv")

        return x



    def test_run_subset(self):
        # computes MPRA score, unactive regions, or other active regions, then concat to df
        x = WeightFeatures.run_subset(sys_args=0, col_names_to_weight_param=0, db_name=0, training_dir_results=0, col_names=0, cell_table=0)
        return

    def test_get_param_weigths(self):
        # executes run_subset, then get_coeff
        x = WeightFeatures.get_param_weights(col_names_to_weight_param=0, db_name=0, motif_info_col_names=0, datafiles_motifs_dir=0,
                      training_dir_results=0, training_dir_Ernst=0, training_dir_Tewhey=0, training_dir_Vockley=0,
                      datafiles_HepG2_geneexpr_dir=0, datafiles_K562_geneexpr_dir=0, datafiles_GM12878_geneexpr_dir=0,
                      datafiles_IMR90_geneexpr_dir=0,
                      cell_table=0)
        return


    def compare_files(self, fileA, fileB):
        with open(fileA, 'r') as a, open(fileB, 'r') as b:
            differ = difflib.Differ()
            for line in differ.compare(a.readlines(), b.readlines()):
                if line[0] in ['-', '+']:
                    return False
        return True

    def test_score_per_pos(self):
        MPRA_tiles_input_file = "InputTestFilesSection3/baseprediction_infile"
        output_file = "InputTestFilesSection3/baseprediction_test_outcome.txt"
        experiment_cell_name = 'HepG2'
        x = WeightFeatures.score_per_pos(MPRA_tiles_input_file, output_file, experiment_cell_name, tile_length=145, region_pos_col=0,
                      region_chr_index=3, cell_name_index=0, region_center_index=4, region_name_sep='_',
                      values_start_index=1, sep='\t')
        assert self.compare_files(output_file, "InputTestFilesSection3/baseprediction_expected_outcome")
        return x


    def test_get_motif_score(self):
        scores_per_bp_input_file = "InputTestFilesSection3/baseprediction_test_outcome.txt"
        motifs_dir = "InputTestFilesSection3/motifs_split_chr"
        motifs_scored_output_file = "InputTestFilesSection3/motif_score_test_outcome.bed"
        x = WeightFeatures.get_motif_scores(scores_per_bp_input_file, motifs_dir, motifs_scored_output_file)
        # TODO: check with expected outcome: looks good, but is there an easier way?
        return
    '''

    def test_get_cell_info_for_motifs(self):
        input_file = "InputTestFilesSection3/motif_score_test_outcome.bed"
        df_output_file = "InputTestFilesSection3/output_combinedP_motifs.df"
        col_names = ['chr', 'motifstart', 'motifend', 'name', 'score', 'pval', 'strand']
        cell_table = "cell_table"
        x = get_cell_info_for_motifs(input_file, db_name="funmotifsdb", cells=['HepG2', 'liver'], df_output_file=df_output_file, col_names = col_names, cell_table=cell_table)
        x.to_csv(df_output_file + '.tsv', sep='\t')
        # TODO: check if input files should create output
        return
'''
    def test_get_Bed(self):
        coordinates_input_file = 0
        input_file = 0
        output_file =0
        x = WeightFeatures.getBed(coordinates_input_file, input_file, output_file)

        return
'''

if __name__ == '__main__':
    unittest.main()
