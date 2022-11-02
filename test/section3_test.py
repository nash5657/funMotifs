"""
Created on Sep 26, 2022

@author: Mark Melzer

Test functions to check section 1 in funMotifsMain.py
"""

import unittest
import sys
import pandas as pd
import statsmodels.api as sm

sys.path.insert(1, '../src/')

import WeightFeatures


class TestSection3(unittest.TestCase):

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


    '''
    def test_get_coeff(self):
        # also inserts dummy variables etc
        # restructures df a bit to make it suitable for sm.Logit function
        x = WeightFeatures.get_coeff(df=0, cols_to_weight=0, outcome_col=0, col_names_to_weight_param=0, dfout_filename=0)
        return

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
    '''

if __name__ == '__main__':
    unittest.main()
