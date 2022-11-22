"""
Created on Sep 26, 2022

@author: Mark Melzer

Test functions to check section 3 in funMotifsMain.py
"""

import unittest
import sys
import os
import pandas as pd
import statsmodels.api as sm
import difflib
from glob import glob

sys.path.insert(1, '../src/')

import WeightFeatures, Utilities
from GetDBData_TrainingSets import get_cell_info_for_motifs, get_cell_info_for_regions


def compare_files(fileA, fileB):
    """
    Function to compare an output file with the expected output to assert its correctness.
    """
    with open(fileA, 'r') as a, open(fileB, 'r') as b:
        differ = difflib.Differ()
        for line in differ.compare(a.readlines(), b.readlines()):
            if line[0] in ['-', '+']:
                print(line)
                return False
    return True


class TestSection3(unittest.TestCase):
    """
    Class to test the different functions in the file WeightFeatures.py
    and the called functions in GetDBData_TrainingSets.py
    Attention: the database used was manually adapted to produce outcome with the given input
    The necessary table (cell_table) is uploaded in csv format as section3_postgres.csv.
    """

    def test_funMotifs_logit(self):
        """
        Test for funMotifslogit()
        To assert its correctness, the same regression model was run on the same data in R. The results are compared and
        should not be differ, by more than 0.001.
        """
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
        """
        Test for get_coeff()
        For assertion the output file is compared with an expected output file.
        """

        # test if the code in get_coeff works (funMotifs_Logit checked before)
        df = pd.read_csv("InputTestFilesSection3/combineddfs.tsv_raw.tsv", sep='\t')
        outcome_col = 'activity_score'
        col_names_to_weight_param = df.columns
        col_names_to_weight_param = list(col_names_to_weight_param[3:])

        dfout_filename = "./InputTestFilesSection3/get_coeff_outfile.tsv"

        x = WeightFeatures.get_coeff(df=df, cols_to_weight=df[col_names_to_weight_param], outcome_col=outcome_col,
                                     col_names_to_weight_param=col_names_to_weight_param,
                                     dfout_filename=dfout_filename)
        assert compare_files(dfout_filename, "InputTestFilesSection3/combineddfs.tsv")

        return x

    def test_score_per_pos(self):

        """
        Test for score_per_pos()
        For assertion the output file is compared with an expected output file.
        """
        MPRA_tiles_input_file = "InputTestFilesSection3/TrainingSets/tile_prom_region/HepG2/input_data/infile.txt"
        output_file = "InputTestFilesSection3/baseprediction_test_outcome.txt"
        experiment_cell_name = 'HepG2'
        x = WeightFeatures.score_per_pos(MPRA_tiles_input_file, output_file, experiment_cell_name, tile_length=145,
                                         region_pos_col=0,
                                         region_chr_index=3, cell_name_index=0, region_center_index=4,
                                         region_name_sep='_',
                                         values_start_index=1, sep='\t')

        assert compare_files(output_file, "InputTestFilesSection3/baseprediction_expected_outcome.txt")
        return

    def test_get_motif_score(self):
        """
        Test for get_motif_score()
        For assertion the output file is compared with an expected output file.
        """
        scores_per_bp_input_file = "InputTestFilesSection3/baseprediction_test_outcome.txt"
        motifs_dir = "InputTestFilesSection3/motifs_split_chr"
        motifs_scored_output_file = "InputTestFilesSection3/motif_score_expected_outcome.bed"
        x = WeightFeatures.get_motif_scores(scores_per_bp_input_file, motifs_dir, motifs_scored_output_file)
        assert compare_files(motifs_scored_output_file, "InputTestFilesSection3/motif_score_expected_outcome.bed")
        return

    def test_get_cell_info_for_motifs(self):
        """
        Test for get_cell_info_for_motifs()
        For assertion the output file is compared with an expected output file.
        """
        input_file = "InputTestFilesSection3/motif_score_expected_outcome.bed"
        df_output_file = "InputTestFilesSection3/output_combinedP_motifs.df"
        col_names = ['chr', 'motifstart', 'motifend', 'name', 'score', 'pval', 'strand']
        cell_table = "cell_table"
        x = get_cell_info_for_motifs(input_file, db_name="funmotifsdb", cells=['liver'], df_output_file=df_output_file,
                                     col_names=col_names, cell_table=cell_table)
        x.to_csv(df_output_file + '.tsv', sep='\t')

        assert compare_files(df_output_file + ".tsv", "InputTestFilesSection3/output_combinedP_motifs.df_expected.tsv")

        return

    def test_get_promoters_of_unactive_genes(self):
        """
        Test for get_promoters_of_unactive_genes()
        For assertion the output file is compared with an expected output file.
        """
        genes_input_file = "InputTestFilesSection3/HepG2_chr10.bed"
        proms_unact_gene_outfile = "InputTestFilesSection3/HepG2_unactive_proms_chr10_expected.bed"
        x = WeightFeatures.get_promoters_of_unactive_genes(genes_input_file=genes_input_file,
                                                           proms_of_unactive_genes_output_file=proms_unact_gene_outfile,
                                                           unactive_expr_thresh=0, gene_expression_index=-1,
                                                           strand_index=5, sep='\t', num_bp_for_prom=1000)
        assert compare_files(proms_unact_gene_outfile, "InputTestFilesSection3/HepG2_unactive_proms_chr10_expected.bed")

        return

    def test_get_cell_info_for_regions(self):
        """
        Test for get_cell_info_for_regions()
        For assertion the output file is compared with an expected output file.
        """
        proms_unactive_genes = "InputTestFilesSection3/HepG2_unactive_proms_chr10_expected.bed"
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
        assert compare_files(prom_df_output_file + '.tsv',
                             "InputTestFilesSection3/HepG2_unactive_proms_motifs.df_expected.tsv")

        return

    def test_get_Bed(self):
        """
        Test for get_Bed()
        For assertion the output file is compared with an expected output file.
        """
        coordinates_input_file = "InputTestFilesSection3/TrainingSets/other_active_region/HepG2/input_data/" \
                                 "infile_coordinates.bed"
        input_file = "InputTestFilesSection3/TrainingSets/other_active_region/HepG2/input_data/infile_expression.txt"
        output_file = "InputTestFilesSection3/Geuv_MGProbes_testbed_outfile.bed"
        x = WeightFeatures.getBed(coordinates_input_file, input_file, output_file)
        expected_outfile = "InputTestFilesSection3/Geuv_MGProbes_expectbed_outfile.bed"
        assert compare_files(output_file, expected_outfile)
        return

    def test_run_subset(self):
        """
        Test for run_subset()
        For assertion the expected outcome data frame is created and compared with the actual outcome.
        """
        training_data_dir = "InputTestFilesSection3/TrainingSets"
        tile_prom_region_path = training_data_dir + '/tile_prom_region'
        prom_unactive_path = training_data_dir + '/prom_unactive'
        other_act_reg_path = training_data_dir + '/other_active_region'
        training_dir_results = training_data_dir + '/training_results/'

        tissue_for_cell_name, cell_name_for_tissue = \
            Utilities.cell_to_tissue_matches('InputTestFilesSection3/TissueCellMatches')
        motif_split_chr = "../datafiles/Motifs/motifs_per_chr"
        tile_prom_region_data = WeightFeatures.get_trainings_data_dirs(tile_prom_region_path, cell_name_for_tissue,
                                                                       tissue_for_cell_name,
                                                                       motif_split_chr)
        prom_unactive_data = WeightFeatures.get_trainings_data_dirs(prom_unactive_path, cell_name_for_tissue,
                                                                    tissue_for_cell_name)
        other_act_reg_data = WeightFeatures.get_trainings_data_dirs(other_act_reg_path, cell_name_for_tissue,
                                                                    tissue_for_cell_name)
        motif_info_col_names = ['chr', 'motifstart', 'motifend', 'name', 'score', 'pval', 'strand']

        col_names_to_weight_param = ['ChromHMM'.lower(), 'DNase__seq'.lower(), 'FANTOM'.lower(),
                                     'NumOtherTFBinding'.lower(), 'RepliDomain'.lower(), 'TFBinding'.lower(),
                                     'TFExpr'.lower(), 'score'.lower(), 'footprints'.lower(), 'cCRE'.lower(),
                                     'IndexDHS'.lower(), 'RegElem'.lower()]
        x = WeightFeatures.run_subset(tile_prom_region_data=tile_prom_region_data,
                                      prom_unactive_data=prom_unactive_data,
                                      other_act_reg_data=other_act_reg_data,
                                      col_names_to_weight_param=col_names_to_weight_param,
                                      db_name='funmotifsdb', db_user_name='mm99', cell_table='cell_table',
                                      col_names=motif_info_col_names, training_dir_results=training_dir_results)
        d = {'score': [11.9703, 11.7222, 14.3636], 'tfexpr': [67.140, 1.902, 0.1597],
             'activity_score': [10.5, '0.0', '1880.82901534565'], 'cellname': ['liver', 'liver', 'liver']}
        df = pd.DataFrame(d)
        z = (x == df)
        for i in z.all():
            assert i is True
        return

    def test_get_param_weights(self):
        """
        Test for get_param_weights()
        For assertion the expected output was computed externally and is compared with the computed output.
        """
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
        tissues_for_cell_name, cell_name_for_tissue =\
            Utilities.cell_to_tissue_matches("InputTestFilesSection3/TissueCellMatches")

        logit_params = WeightFeatures.get_param_weights(training_data_dir, col_names_to_weight_param, db_name,
                                                        motif_info_col_names, cell_table, db_user_name,
                                                        cell_name_for_tissue, tissues_for_cell_name,
                                                        motif_split_chr=datafiles_motifs_dir)

        assert logit_params.params[0] - 0.0576878 < 0.001
        assert logit_params.params[1] - 0.2165892 < 0.001

        return

    def test_get_trainings_data_dir(self):
        """
        Test for get_trainings_data_dir()
        For assertion the expected output is compared with the actual output.
        """
        tissue_for_cell_name, cell_name_for_tissue = \
            Utilities.cell_to_tissue_matches("InputTestFilesSection3/TissueCellMatches")
        motif_split_chr = "../datafiles/Motifs/motifs_per_chr"
        tile_prom_region_path = "InputTestFilesSection3/TrainingSets/tile_prom_region"

        x = WeightFeatures.get_trainings_data_dirs(tile_prom_region_path, cell_name_for_tissue, tissue_for_cell_name,
                                                   motif_split_chr)

        assert x == [['InputTestFilesSection3/TrainingSets/tile_prom_region/HepG2/input_data/infile.txt',
                      'InputTestFilesSection3/TrainingSets/tile_prom_region/HepG2/output_data', 'HepG2',
                      ['HepG2', 'Liver'], '../datafiles/Motifs/motifs_per_chr']]
        return

    def test_tile_prom_region(self):
        """
        Test for tile_prom_region()
        For assertion the expected output is compared with the actual output.
        """
        col_names = ['chr', 'motifstart', 'motifend', 'name', 'score', 'pval', 'strand']
        cell_table = "cell_table"
        path = "InputTestFilesSection3/TrainingSets/tile_prom_region"
        tissue_for_cell_name, cell_name_for_tissue = \
            Utilities.cell_to_tissue_matches("InputTestFilesSection3/TissueCellMatches")
        motif_split_chr = "../datafiles/Motifs/motifs_per_chr"
        data = WeightFeatures.get_trainings_data_dirs(path, cell_name_for_tissue, tissue_for_cell_name,
                                                      motif_split_chr)[0]
        MPRA_score_per_pos_outfile = data[1] + '/MPRA_score_per_pos_outfile.txt'
        motifs_scored_output_file = data[1] + '/motifs_scored_output_file.bed'
        df_output_file = data[1] + '/motifs_scored_output_dataframe.df'
        x = WeightFeatures.tile_prom_regions(MPRA_infile=data[0],
                                             MPRA_score_per_pos_outfile=MPRA_score_per_pos_outfile,
                                             experiment_cell_name=data[2], cells_to_extract_info_from=data[3],
                                             motifs_dir=data[4],
                                             motifs_scored_output_file=motifs_scored_output_file,
                                             df_output_file=df_output_file, db_name='funmotifsdb',
                                             db_user_name='mm99', col_names=col_names,
                                             cell_table=cell_table)
        d = {'chr': 10, 'motifstart': 21798914, 'motifend': 21798927, 'name': 'PLAG1_MA0163.1', 'score': 11.9703,
             'pval': 0.0000651, 'strand': '-', 'liver___tfexpr': 67.14, 'Activity_Score': 10.5}
        df = pd.DataFrame(d, index=[0])
        z = (x == df)
        for i in z.all():
            assert i is True
        return

    def test_prom_unactive(self):
        """
        Test for prom_unactive()
        For assertion the expected output is compared with the actual output.
        """
        path = "InputTestFilesSection3/TrainingSets/prom_unactive"
        tissue_for_cell_name, cell_name_for_tissue = \
            Utilities.cell_to_tissue_matches("InputTestFilesSection3/TissueCellMatches")
        motif_split_chr = "../datafiles/Motifs/motifs_per_chr"
        data = WeightFeatures.get_trainings_data_dirs(path, cell_name_for_tissue, tissue_for_cell_name,
                                                      motif_split_chr)[0]
        proms_unact_genes_outfile = data[1] + '/' + data[0].split('/')[-1].split('.')[0] + '_unactive_proms.bed'
        prom_df_outfile = data[1] + '/' + data[0].split('/')[-1].split('.')[0] + '_unactive_proms.df'
        col_names = ['chr', 'motifstart', 'motifend', 'name', 'score', 'pval', 'strand']

        x = WeightFeatures.prom_unactive(geneexp_infile=data[0],
                                         proms_unact_genes_outfile=proms_unact_genes_outfile,
                                         prom_df_outfile=prom_df_outfile, cells_to_extract_info_from=data[2],
                                         db_name='funmotifsdb', db_user_name='mm99', col_names=col_names)

        d = {'chr': 10, 'motifstart': 50116528, 'motifend': 50117528, 'name': 'PPARG_MA0066.1', 'score': 11.7222,
             'pval': 0.0000104, 'strand': '-', 'Activity_Score': '0.0'}
        df = pd.DataFrame(d, index=[0])
        z = (x == df)
        for i in z.all():
            assert i is True
        return

    def test_other_active_region(self):
        """
        Test for other_active_region()
        For assertion the expected output is compared with the actual output.
        """
        path = "InputTestFilesSection3/TrainingSets/other_active_region"
        tissue_for_cell_name, cell_name_for_tissue = \
            Utilities.cell_to_tissue_matches("InputTestFilesSection3/TissueCellMatches")
        motif_split_chr = "../datafiles/Motifs/motifs_per_chr"
        data = WeightFeatures.get_trainings_data_dirs(path, cell_name_for_tissue, tissue_for_cell_name,
                                                      motif_split_chr)[0]
        other_act_reg_file = data[1] + '/output.bed'
        regions_df_outfile = data[1] + '/output_motifs.df'
        col_names = ['chr', 'motifstart', 'motifend', 'name', 'score', 'pval', 'strand']

        x = WeightFeatures.other_active_region(coordinates_infile=data[0], infile=data[2],
                                               other_act_reg_file=other_act_reg_file,
                                               cells_to_extract_info_from=data[3],
                                               regions_df_outfile=regions_df_outfile, db_name='funmotifsdb',
                                               db_user_name='mm99', col_names=col_names)

        d = {'chr': 10, 'motifstart': 101769701, 'motifend': 101769711, 'name': 'FOSL1_MA0477.1', 'score': 14.3636,
             'pval': 0.000013, 'strand': '+', 'liver___tfexpr': 0.1597, 'Activity_Score': '1880.82901534565'}
        df = pd.DataFrame(d, index=[0])
        z = (x == df)
        for i in z.all():
            assert i is True

        return


if __name__ == '__main__':
    prev = glob("InputTestFilesSection3/*")
    unittest.main(exit=False)
    rm = glob("InputTestFilesSection3/TrainingSets/*/*/output_data/*")
    for i in glob("InputTestFilesSection3/TrainingSets/training_results/*"):
        rm.append(i)
    after = glob("InputTestFilesSection3/*")
    for i in after:
        if i not in prev:
            rm.append(i)
    for i in rm:
        os.remove(i)
