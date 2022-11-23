"""
Created on Dec 9, 2016

@author: husensofteng
@contributor: markmelzer
"""

import os
from glob import glob
import numpy as np
import pandas as pd
import statsmodels.api as sm
from pybedtools import BedTool

from GetDBData_TrainingSets import get_cell_info_for_motifs, get_cell_info_for_regions


def score_per_pos(MPRA_tiles_input_file, output_file, experiment_cell_name, tile_length=145, region_pos_col=0,
                  region_chr_index=3, cell_name_index=0, region_center_index=4, region_name_sep='_',
                  values_start_index=1, sep='\t'):
    with open(MPRA_tiles_input_file, 'r') as infile, open(output_file, 'w') as outfile:
        line = infile.readline()
        while line:
            out_lines = []
            sline = line.split()
            region_chr = sline[region_pos_col].split(region_name_sep)[region_chr_index]
            region_center = sline[region_pos_col].split(region_name_sep)[region_center_index]
            values = sline[values_start_index::]
            curr_pos = tile_length + 2
            before_center = True
            for v in values:
                if curr_pos == 0:
                    before_center = False
                if before_center:  # set the start and end for each region
                    out_lines.append([region_chr, str(int(region_center) - (tile_length - curr_pos)),
                                      str(int(region_center) - (tile_length - curr_pos)), experiment_cell_name, v,
                                      str(curr_pos), sline[region_pos_col]])
                    curr_pos -= 1
                else:
                    out_lines.append(
                        [region_chr, str(int(region_center) + curr_pos), str(int(region_center) + curr_pos),
                         experiment_cell_name, v, str(curr_pos), sline[region_pos_col]])
                    curr_pos += 1
            line = infile.readline()
            for out_line in out_lines:
                outfile.write(sep.join(out_line) + '\n')
    return


def get_motif_scores(scores_per_bp_input_file, motifs_dir, motifs_scored_output_file, max_motif_score_per_region=False):
    regions = {}
    scores_per_bp_input_file_obj = BedTool(scores_per_bp_input_file)
    for motif_file in os.listdir(motifs_dir):
        motif_file_obj = BedTool(motifs_dir + '/' + motif_file)
        scores_per_bp_input_file_results_tmp_inter = scores_per_bp_input_file + '_' + motif_file
        # TODO: next command messes up the nummeration: compare scores_per_bp_input_file_results_tmp_inter
        #  and motifs_scored_output_file; https://daler.github.io/pybedtools/intersections.html
        motif_file_obj.intersect(scores_per_bp_input_file_obj, wo=True).saveas(
            scores_per_bp_input_file_results_tmp_inter)

        if os.stat(scores_per_bp_input_file_results_tmp_inter).st_size != 0:

            scores_per_bp_input_file_results_tmp_inter_sorted = scores_per_bp_input_file_results_tmp_inter + '_sorted'
            os.system("""sort -V -k1,1 -k2n,2 -k3n,3 -k4,4 -k5,5 -k6,6 -k7,7 {} > {}""".format(
                scores_per_bp_input_file_results_tmp_inter, scores_per_bp_input_file_results_tmp_inter_sorted))

            results = BedTool(scores_per_bp_input_file_results_tmp_inter_sorted).groupby(g=[1, 2, 3, 4, 5, 6, 7],
                                                                                         c=[11, 12, 13, 12, 14],
                                                                                         o=['distinct', 'mean', 'min',
                                                                                            'collapse', 'distinct'])
            scores_per_bp_input_file_results_tmp = scores_per_bp_input_file_results_tmp_inter_sorted + '_tmp'
            with open(results.fn, 'r') as r, open(scores_per_bp_input_file_results_tmp, 'w') as results_outfile:
                line = r.readline()
                while line:
                    sl = line.strip().split('\t')
                    results_outfile.write(
                        '\t'.join(sl[0:4]) + '\t' + 'S' + sl[4] + 'P' + sl[5] + '\t' + '\t'.join(sl[6:12]) + '\n')
                    line = r.readline()

            with open(scores_per_bp_input_file_results_tmp, 'r') as r:
                line = r.readline()
                while line:
                    sl = line.strip().split('\t')
                    if len(sl[9].split(',')) == (
                            (int(sl[2]) - int(sl[1])) + 1):  # only keep if all bases of the motif were in the region
                        if sl[10] not in list(regions.keys()):
                            regions[sl[10]] = [sl[0:10]]
                        else:
                            if max_motif_score_per_region:
                                if float(sl[7]) >= float(regions[sl[10]][0][7]):
                                    regions[sl[10]] = [sl[0:10]]
                            else:
                                regions[sl[10]].append(sl[0:10])
                    line = r.readline()

    with open(motifs_scored_output_file, 'w') as motifs_scored_outfile:
        for region in list(regions.keys()):
            for motif_region in regions[region]:
                motif_region[7] = "{:.3f}".format(float(motif_region[7]))
                motifs_scored_outfile.write('\t'.join(motif_region) + '\t' + region + '\n')

        # TODO: check try except structure belo: check try except structure below
        try:
            os.remove(scores_per_bp_input_file_results_tmp)
        except:
            pass
    return regions  # the highest scored motif for each tile region


def get_promoters_of_unactive_genes(genes_input_file, proms_of_unactive_genes_output_file, unactive_expr_thresh=0,
                                    gene_expression_index=-1, strand_index=4, sep='\t', num_bp_for_prom=1000):
    # if os.path.exists(proms_of_unactive_genes_output_file):
    #    return proms_of_unactive_genes_output_file
    with open(genes_input_file, 'r') as genes_infile, open(proms_of_unactive_genes_output_file, 'w') as outfile:
        l = genes_infile.readline()
        while l:
            sl = l.strip().split(sep)
            if float(sl[gene_expression_index]) <= unactive_expr_thresh and sl[8] == 'protein_coding':
                if sl[strand_index] == "-":
                    outfile.write(sl[0] + sep + sl[2] + sep + str(int(sl[2]) + num_bp_for_prom) + sep + sep.join(
                        sl[3:-1]) + sep + '0.0' + '\n')
                else:
                    outfile.write(sl[0] + sep + str(int(sl[1]) - num_bp_for_prom) + sep + sl[1] + sep + sep.join(
                        sl[3:-1]) + sep + '0.0' + '\n')
            l = genes_infile.readline()
    return proms_of_unactive_genes_output_file


def getBed(coordinates_input_file, input_file, output_file):
    with open(coordinates_input_file, 'r') as oligos:

        oligos_dict = {}
        for o in oligos:
            so = o.strip().split('\t')
            if so[3] not in list(oligos_dict.keys()):
                oligos_dict[so[3]] = [so[0], so[1], so[2]]

    with open(input_file, 'r') as lines, open(output_file, 'w') as out_file:

        for o in lines:
            so = o.strip().split(' ')
            # TODO: check why so[1] might not be in oligos_dict.keys()
            if so[3] == 'ref' and (float(so[8]) >= 2 or float(so[13]) >= 2) and so[1] in oligos_dict.keys():
                out_file.write(
                    '\t'.join(oligos_dict[so[1]]) + '\t' + so[0] + '\t' + str(max(float(so[5]), float(so[10]))) + '\n')

    return output_file


def tile_prom_regions(MPRA_infile, MPRA_score_per_pos_outfile, experiment_cell_name, cells_to_extract_info_from,
                      motifs_dir, motifs_scored_output_file, df_output_file, db_name, col_names, cell_table,
                      db_user_name, max_motif_score_per_region=False, sep='\t'):
    """
        Get MPRA (Activity) score for the motifs
    """
    #if not os.path.exists(MPRA_score_per_pos_outfile):
    score_per_pos(MPRA_tiles_input_file=MPRA_infile, output_file=MPRA_score_per_pos_outfile,
                  experiment_cell_name=experiment_cell_name)

    #if not os.path.exists(motifs_scored_output_file):
    get_motif_scores(scores_per_bp_input_file=MPRA_score_per_pos_outfile,
                     motifs_dir=motifs_dir, motifs_scored_output_file=motifs_scored_output_file,
                     max_motif_score_per_region=max_motif_score_per_region)

    df_results = get_cell_info_for_motifs(motifs_scored_output_file, db_name=db_name,
                                          cells=cells_to_extract_info_from, df_output_file=df_output_file,
                                          col_names=col_names, cell_table=cell_table, db_user_name=db_user_name)
    df_results.to_csv(df_output_file + '.tsv', sep=sep)

    return df_results


def prom_unactive(geneexp_infile, proms_unact_genes_outfile, prom_df_outfile, cells_to_extract_info_from, db_name,
                  db_user_name, col_names, cols_indices_to_report_from_file=[9], assays=['all'], unactive_expr_thresh=0,
                  gene_expression_index=-1, strand_index=5, sep='\t', num_bp_for_prom=1000, cell_table='cell_table'):
    """
        Get those promoters that are inactive, i.e. promoters to unexpressed genes
    """
    proms_unactive_genes = get_promoters_of_unactive_genes(genes_input_file=geneexp_infile,
                                                           proms_of_unactive_genes_output_file=proms_unact_genes_outfile,
                                                           unactive_expr_thresh=unactive_expr_thresh,
                                                           gene_expression_index=gene_expression_index,
                                                           strand_index=strand_index, sep=sep,
                                                           num_bp_for_prom=num_bp_for_prom)
    prom_results_df = get_cell_info_for_regions(proms_unactive_genes, db_user_name=db_user_name, db_name=db_name,
                                                cells=cells_to_extract_info_from, assays=assays, cell_table=cell_table,
                                                df_output_file=prom_df_outfile, col_names=col_names,
                                                cols_indices_to_report_from_file=cols_indices_to_report_from_file)
    prom_results_df.to_csv(prom_df_outfile + '.tsv', sep=sep)
    return prom_results_df


def other_active_region(coordinates_infile, infile, other_act_reg_file, cells_to_extract_info_from, regions_df_outfile,
                        db_name, db_user_name, col_names, assays=['all'], sep='\t', cell_table='cell_table',
                        cols_indices_to_report_from_file=[4],
                        cols_names_to_report_from_file=['Activity_Score'], region_name_index=3,
                        region_strand_index=None, region_score_index=4, motif_score_index=4,
                        max_number_motifs_to_report=1, min_dist_from_region_start=False,
                        min_dist_from_region_center=True, max_motif_score=False, max_region_score=False):
    """
        Find other active regions in the cell
    """
    other_act_reg_file = getBed(coordinates_input_file=coordinates_infile, input_file=infile,
                                    output_file=other_act_reg_file)


    other_active_regions_df = get_cell_info_for_regions(
        other_act_reg_file, db_user_name=db_user_name, db_name=db_name, cells=cells_to_extract_info_from, assays=assays,
        sep=sep, cols_indices_to_report_from_file=cols_indices_to_report_from_file,
        cols_names_to_report_from_file=cols_names_to_report_from_file, df_output_file=regions_df_outfile,
        region_name_index=region_name_index, region_strand_index=region_strand_index,
        region_score_index=region_score_index, motif_score_index=motif_score_index,
        max_number_motifs_to_report=max_number_motifs_to_report, min_dist_from_region_start=min_dist_from_region_start,
        min_dist_from_region_center=min_dist_from_region_center, max_motif_score=max_motif_score,
        max_region_score=max_region_score, col_names=col_names, cell_table=cell_table)
    other_active_regions_df.to_csv(regions_df_outfile + '.tsv', sep=sep)
    return other_active_regions_df


def run_subset(tile_prom_region_data: list, prom_unactive_data: list, other_act_reg_data: list,
               col_names_to_weight_param, db_name, training_dir_results, col_names, cell_table, db_user_name):
    """
        Gather the data to perform the regression on it
    """

    for data in tile_prom_region_data:
        MPRA_score_per_pos_outfile = data[1] + '/MPRA_score_per_pos_outfile.txt'
        motifs_scored_output_file = data[1] + '/motifs_scored_output_file.bed'
        df_output_file = data[1] + '/motifs_scored_output_dataframe.df'

        df1 = tile_prom_regions(MPRA_infile=data[0],
                                MPRA_score_per_pos_outfile=MPRA_score_per_pos_outfile,
                                experiment_cell_name=data[2], cells_to_extract_info_from=data[3],
                                motifs_dir=data[4],
                                motifs_scored_output_file=motifs_scored_output_file,
                                df_output_file=df_output_file, db_name=db_name,
                                db_user_name=db_user_name, col_names=col_names,
                                cell_table=cell_table)

    for data in prom_unactive_data:
        proms_unact_genes_outfile = data[1] + '/' + data[0].split('/')[-1].split('.')[0] + '_unactive_proms.bed'
        prom_df_outfile = data[1] + '/' + data[0].split('/')[-1].split('.')[0] + '_unactive_proms.df'

        df2 = prom_unactive(geneexp_infile=data[0],
                            proms_unact_genes_outfile=proms_unact_genes_outfile,
                            prom_df_outfile=prom_df_outfile, cells_to_extract_info_from=data[2],
                            db_name=db_name, db_user_name=db_user_name, col_names=col_names, cell_table=cell_table)

    for data in other_act_reg_data:
        other_act_reg_file = data[1] + '/output.bed'
        regions_df_outfile = data[1] + '/output_motifs.df'

        df3 = other_active_region(coordinates_infile=data[0], infile=data[2],
                                  other_act_reg_file=other_act_reg_file,
                                  cells_to_extract_info_from=data[3],
                                  regions_df_outfile=regions_df_outfile, db_name=db_name,
                                  db_user_name=db_user_name, col_names=col_names, cell_table=cell_table)
    # TODO: check difference and use of prom_unactive and other_active_region
    # combine data frames
    combined_results = pd.concat([df1, df2, df3])
    # create data frame for return values
    dfout = pd.DataFrame()
    dfout_filename = "dfout_file"
    cell_name = ""

    # get columns that should be weighted
    if len(combined_results) > 0:
        outcome_col = 'Activity_Score'
        cols_to_weight = []
        for c in combined_results.columns:
            if c.split('___')[-1] in col_names_to_weight_param:
                cols_to_weight.append(c)

        # TODO: why break? what happens when more than 1 cell in data or is this not possible?
        # does it only determine cell name? then there should be a better way
        for k in cols_to_weight:
            if '___' in k:
                dfout_filename = k.split('___')[0] + "combineddfs.tsv"
                cell_name = k.split('___')[0]
                break

        # add columns that shall be weighted into data frame
        dfout = combined_results.loc[:, cols_to_weight]
        dfout[outcome_col] = combined_results.loc[:, outcome_col]

        # add cell names to data frame
        cell_name_series = pd.Series([cell_name for x in range(0, len(dfout))], name='CellName')
        dfout['CellName'] = cell_name_series
        dfout_cols = {}

        for c in dfout.columns:
            if '___' in c:
                dfout_cols[c] = c.split('___')[-1].lower()
            else:
                dfout_cols[c] = c.lower()
                dfout_cols[c.encode('ascii', 'ignore')] = c.encode('ascii', 'ignore').lower()

        # format data frame and convert to csv
        if not os.path.exists(training_dir_results):
            os.makedirs(training_dir_results)
        dfout.rename(columns=dfout_cols, inplace=True)
        dfout.to_csv(training_dir_results + dfout_filename, sep='\t')

    dfout['activity_score'] = dfout['activity_score'].astype(float)

    return dfout.reset_index(drop=True)


def make_binary(x, f=0):
    if type(x) is str:
        return x
    if x > f:
        return 1
    else:
        return 0


def make_abs(x):
    if type(x) is str:
        return x
    return abs(x)


def inf_to_zero(x):
    if abs(x) == np.inf:
        return 0.0
    else:
        return x


def combine_lables(x):
    if 'Tss' in x or 'Tx' in x or 'BivFlnk' in x:
        return 'TSS'
    elif 'Enh' in x:
        return 'Enh'
    elif 'Repr' in x:
        return 'Repr'
    elif 'Quies' in x or 'Rpts' in x or 'Het' in x:
        return 'Quies'
    else:
        return x


def get_coeff(df, cols_to_weight, outcome_col, col_names_to_weight_param, dfout_filename):
    """
    Compute coefficients for the logistic regression model
    """
    # create new data frame with output column
    new_df = pd.DataFrame()
    # TODO: is the order of the applys below right? (threshold for binary vs absolute)
    new_df[outcome_col] = df[outcome_col].astype(float).apply(make_abs).apply(make_binary, args=(0,))
    # TODO: isn't this unnecessary???
    new_df[outcome_col] = new_df[outcome_col]

    # prepare values for parameter calculation (make binary etc.)
    for c in cols_to_weight:
        if c not in col_names_to_weight_param:
            continue
        if c == 'DNase__seq'.lower() or c == 'TFBinding'.lower() or c == 'FANTOM'.lower() or c == 'footprints'.lower():
            new_df[c] = df[c].apply(make_binary, args=(0,))
        elif c == 'NumOtherTFBinding'.lower():
            new_df[c] = df[c]
        elif c == 'RepliDomain'.lower() or c == 'CellName'.lower() or c == 'name' or c == 'IndexDHS'.lower() or \
                c == 'cCRE'.lower() or c == 'RegElem'.lower():
            df_dummies = pd.get_dummies(df[c])
            new_df = pd.concat([new_df, df_dummies], axis=1)
        elif c == 'ChromHMM'.lower():
            df[c] = df[c].apply(combine_lables)
            df_dummies = pd.get_dummies(df[c])
            new_df = pd.concat([new_df, df_dummies], axis=1)
        elif c == 'score'.lower():
            new_df[c] = df[c]
        elif c == 'TFExpr'.lower():
            new_df[c] = df[c].apply(float).apply(np.log10).apply(inf_to_zero)

    new_cols_to_weight = [c for c in new_df.columns]

    # new_cols_to_weight = [c.encode('ascii','ignore') for c in new_df.columns]

    while 'NO' in new_cols_to_weight:  # remove NO == no overlap labels
        del new_cols_to_weight[new_cols_to_weight.index('NO')]

    df.to_csv(dfout_filename + '_raw.tsv', sep='\t')
    new_df.to_csv(dfout_filename, sep='\t')

    del new_cols_to_weight[new_cols_to_weight.index(outcome_col)]

    # compute parameters & return results
    return dfout_filename, funMotifs_logit(new_df[outcome_col], new_df[new_cols_to_weight])


def funMotifs_logit(outcome_col, col_weight, method='bfgs', maxiter=10000, full_output=True):
    # TODO: this deletes all columns containing NaN's, but figure out why there are NaNs in the first place
    col_weight = col_weight.dropna(axis=1)
    # TODO: remove prints above and check how to incoporate the .astype(float) below best
    model = sm.Logit(outcome_col.astype(float), col_weight.astype(float)).fit(method=method, maxiter=maxiter,
                                                                              full_output=full_output)

    return model


def get_trainings_data_dirs(path, cell_name_for_tissue, tissue_for_cell_name, motif_split_chr=''):
    # TODO: write unit test
    # only for tile_prom_regions so far
    cell_list = glob(path + '/*')
    data_lists = []
    for cell in cell_list:
        try:
            data_list = [None] * 5
            # get input files as first input for all three cases
            if path.__contains__('tile'):
                data_list[0] = cell + '/input_data/infile.txt'
            elif path.__contains__('unact'):
                data_list[0] = cell + '/input_data/' + cell.split('/')[-1] + '.bed'
            elif path.__contains__('other'):
                data_list[0] = cell + '/input_data/infile_coordinates.bed'
            else:
                continue

            # create or set output directory for all three cases
            data_list[1] = cell + '/output_data'
            if not os.path.exists(data_list[1]):
                os.makedirs(data_list[1])

            # individual inputs for all three cases
            if path.__contains__('unact'):
                try:
                    data_list[2] = tissue_for_cell_name[cell_name_for_tissue[cell.split('/')[-1]]]
                except:
                    data_list[2] = [cell.split('/')[-1]]
                    pass
                # done for unactive promoters
                data_lists.append(data_list)
                continue
            elif path.__contains__('other'):
                data_list[2] = cell + '/input_data/infile_expression.txt'
            elif path.__contains__('tile'):
                data_list[2] = cell.split('/')[-1]
            else:
                continue

            # get representative cell name for all cells of tissue

            try:
                data_list[3] = cell_name_for_tissue[tissue_for_cell_name[cell.split('/')[-1]]]
            except:
                data_list[3] = [cell.split('/')[-1]]
                pass

            # add motif path for tile promoter region
            if path.__contains__('tile'):
                data_list[4] = motif_split_chr
            data_lists.append(data_list)
        except Exception as e:
            print("Could not determine trainings data for " + cell.split('/')[-1])
            print(e)
            pass

    return data_lists


def get_param_weights(training_data_dir, col_names_to_weight_param, db_name, motif_info_col_names, cell_table,
                      db_user_name, cell_name_for_tissue, tissue_for_cell_name,
                      motif_split_chr="../datafiles/Motifs/motifs_per_chr"):
    # TODO: introduce force-overwrite to this file

    tile_prom_region_path = training_data_dir + '/tile_prom_region'
    prom_unactive_path = training_data_dir + '/prom_unactive'
    other_act_reg_path = training_data_dir + '/other_active_region'
    training_dir_results = training_data_dir + '/training_results/'

    tile_prom_region_data = get_trainings_data_dirs(tile_prom_region_path, cell_name_for_tissue, tissue_for_cell_name,
                                                    motif_split_chr)
    prom_unactive_data = get_trainings_data_dirs(prom_unactive_path, cell_name_for_tissue, tissue_for_cell_name)
    other_act_reg_data = get_trainings_data_dirs(other_act_reg_path, cell_name_for_tissue, tissue_for_cell_name)

    df = run_subset(tile_prom_region_data=tile_prom_region_data, prom_unactive_data=prom_unactive_data,
                    other_act_reg_data=other_act_reg_data, col_names_to_weight_param=col_names_to_weight_param,
                    db_name=db_name, db_user_name=db_user_name, cell_table=cell_table,
                    col_names=motif_info_col_names, training_dir_results=training_dir_results)

    # combine training files
    # print(len(df))
    # combined_results = df[0].append(df[1::], ignore_index=True)
    combined_results = df

    # TODO: check why this is done in this way
    reported_col_names = combined_results.columns
    outcome_col = 'Activity_Score'.lower()
    cols_to_weight = [c for c in reported_col_names]
    # cols_to_weight = [c.encode('ascii','ignore') for c in reported_col_names]
    del cols_to_weight[cols_to_weight.index(outcome_col)]

    # prepare output file
    dfout_filename = "{0}/combineddfs.tsv".format(training_dir_results)
    combined_results.to_csv('{0}/combined_results.tsv'.format(training_dir_results), sep='\t')

    # compute coefficients of log model
    dfout_filename, logit_params = get_coeff(combined_results, cols_to_weight, outcome_col, col_names_to_weight_param,
                                             dfout_filename=dfout_filename)

    return logit_params
