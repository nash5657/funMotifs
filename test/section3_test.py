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

import WeightFeatures, Utilities
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
        MPRA_tiles_input_file = "InputTestFilesSection3/infile.txt"
        output_file = "InputTestFilesSection3/baseprediction_test_outcome.txt"
        experiment_cell_name = 'HepG2'
        x = WeightFeatures.score_per_pos(MPRA_tiles_input_file, output_file, experiment_cell_name, tile_length=145, region_pos_col=0,
                      region_chr_index=3, cell_name_index=0, region_center_index=4, region_name_sep='_',
                      values_start_index=1, sep='\t')
        #assert self.compare_files(output_file, "InputTestFilesSection3/baseprediction_expected_outcome")
        return



    def test_get_motif_score(self):
        scores_per_bp_input_file = "InputTestFilesSection3/baseprediction_test_outcome.txt"
        motifs_dir = "InputTestFilesSection3/motifs_split_chr"
        motifs_scored_output_file = "InputTestFilesSection3/motif_score_test_outcome.bed"
        x = WeightFeatures.get_motif_scores(scores_per_bp_input_file, motifs_dir, motifs_scored_output_file)
        # TODO: check with expected outcome: looks good, but is there an easier way?
        return



    def test_get_cell_info_for_motifs(self):
        input_file = "InputTestFilesSection3/motif_score_test_outcome.bed"
        df_output_file = "InputTestFilesSection3/output_combinedP_motifs.df"
        col_names = ['chr', 'motifstart', 'motifend', 'name', 'score', 'pval', 'strand']
        cell_table = "cell_table"
        x = get_cell_info_for_motifs(input_file, db_name="funmotifsdb", cells=['brain'], df_output_file=df_output_file, col_names = col_names, cell_table=cell_table)
        x.to_csv(df_output_file + '.tsv', sep='\t')
        # TODO: check if input files should create output
        return



    def test_get_promoters_of_unactive_genes(self):
        genes_input_file = "InputTestFilesSection3/HepG2.bed"
        proms_unact_gene_outfile = "InputTestFilesSection3/HepG2_unactive_proms.bed"
        x = WeightFeatures.get_promoters_of_unactive_genes(genes_input_file=genes_input_file,
                                                           proms_of_unactive_genes_output_file=proms_unact_gene_outfile,
                                                           unactive_expr_thresh=0, gene_expression_index=-1,
                                                           strand_index=5, sep='\t', num_bp_for_prom=1000)

        return



    def test_get_cell_info_for_regions(self):
        proms_unactive_genes = "InputTestFilesSection3/HepG2_unactive_proms_chr10.bed"  # output of get_proms_of_unactive_genes
        db_user_name = "mm99"
        db_name = "funmotifsdb"
        cells_to_extract_info_from_for_prom_unactive = ['all']
        prom_df_output_file = "InputTestFilesSection3/HepG2_unactive_proms_motifs.df"
        col_names = ['chr', 'motifstart', 'motifend', 'name', 'score', 'pval', 'strand']

        x = WeightFeatures.get_cell_info_for_regions(proms_unactive_genes, db_user_name=db_user_name, db_name=db_name,
                                                     cells=cells_to_extract_info_from_for_prom_unactive, assays=['all'],
                                                     df_output_file=prom_df_output_file, col_names=col_names,
                                                     cols_indices_to_report_from_file=[8])
        x.to_csv(prom_df_output_file + '.tsv', sep='\t')
        # TODO: compare output

        return


    def test_get_Bed(self):
        coordinates_input_file = "InputTestFilesSection3/infile_coordinates.bed"
        input_file = "InputTestFilesSection3/infile_expression.txt"
        output_file = "InputTestFilesSection3/Geuv_MGProbes_testbed_outfile.bed"
        x = WeightFeatures.getBed(coordinates_input_file, input_file, output_file)
        # TODO: compare output
        return


    def test_run_subset(self):
        # TODO: worked before but changed structure
        # computes MPRA score, unactive regions, or other active regions, then concat to df
        sys_args = ["InputTestFilesSection3/infile.txt",
                    "InputTestFilesSection3/baseprediction_test_outcome",
                    "HepG2", "HepG2,liver", "InputTestFilesSection3/motifs_split_chr/",
                    "InputTestFilesSection3/motif_score_test_outcome.bed",
                    "InputTestFilesSection3/motif_score_test_outcome.df",
                    "InputTestFilesSection3/HepG2_chr10.bed",
                    "InputTestFilesSection3/HepG2_unactive_proms_chr10.bed",
                    "InputTestFilesSection3/HepG2_unactive_proms_chr10.df",
                    "HepG2,liver", "InputTestFilesSection3/infile_coordinates.bed",
                    "InputTestFilesSection3/infile_expression.txt",
                    "InputTestFilesSection3/Geuv_MGProbes_testbed_outfile.bed",
                    "HepG2,liver", "InputTestFilesSection3/HepG2Motifs_Tewhy",
                    "tile_prom_regions,prom_unactive,other_active_regions"]
        col_names_to_weight_param = ['ChromHMM'.lower(), 'DNase__seq'.lower(), 'FANTOM'.lower(),
                                     'NumOtherTFBinding'.lower(), 'RepliDomain'.lower(), 'TFBinding'.lower(),
                                     'TFExpr'.lower(), 'score'.lower(), 'footprints'.lower(), 'cCRE'.lower(),
                                     'IndexDHS'.lower(), 'RegElem'.lower()]
        motif_info_col_names = ['chr', 'motifstart', 'motifend', 'name', 'score', 'pval', 'strand']

        x, y = WeightFeatures.run_subset(sys_args=sys_args, col_names_to_weight_param=col_names_to_weight_param,
                                         db_name='funmotifsdb', db_user_name='mm99',
                                         training_dir_results="InputTestFilesSection3/training_dir_result",
                                         col_names=motif_info_col_names, cell_table='cell_table')
        print(x, y)
        return
    '''

    def test_get_param_weights(self):
        cell_table = 'cell_table'
        datafiles_motifs_dir = 'InputTestFilesSection3/motifs_split_chr'

        training_data_dir = 'InputTestFilesSection3/TrainingSets'

        motif_info_col_names = ['chr', 'motifstart', 'motifend', 'name', 'score', 'pval', 'strand']

        col_names_to_weight_param = ['ChromHMM'.lower(), 'DNase__seq'.lower(), 'FANTOM'.lower(),
                                     'NumOtherTFBinding'.lower(), 'RepliDomain'.lower(), 'TFBinding'.lower(),
                                     'TFExpr'.lower(), 'score'.lower(), 'footprints'.lower(), 'cCRE'.lower(),
                                     'IndexDHS'.lower(), 'RegElem'.lower()]
        db_name = 'funmotifsdb'
        db_user_name = 'mm99'
        tissues_for_cell_name, cell_name_for_tissue = Utilities.cell_to_tissue_matches("../conf/TissueCellMatches")

        logit_params = WeightFeatures.get_param_weights(training_data_dir, col_names_to_weight_param, db_name,
                                                        motif_info_col_names, cell_table, db_user_name,
                                                        cell_name_for_tissue, tissues_for_cell_name,
                                                        motif_split_chr=datafiles_motifs_dir)

        print("Result:")
        print(logit_params.params)

        return

if __name__ == '__main__':
    unittest.main()
