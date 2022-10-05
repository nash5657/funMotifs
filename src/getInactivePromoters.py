"""
Created on 04 Oct 2022

@author: Mark Melzer
"""

'''
Extract promoters that are called inactive (see in paper) as training data for the logistic regression model
'''

import sys, os
import pandas as pd
from DBUtilities import open_connection, close_connection, get_col_names_from_table


def updateColNames(col_names=[]):
    for i in range(0, len(col_names)):
        col_names[i] = ('_'.join((col_names[i].replace('(', '').replace(')', '').replace('-', '__')).split())).lower()
    return col_names


def get_cond_stmt(conds=[]):
    if len(conds) > 0:
        return ' where {} '.format(' and '.join(conds))
    else:
        return ''


def get_limit_smt(number_rows_select='all'):
    limit_number_rows_select_stmt = ""
    if number_rows_select != "all":
        if number_rows_select > 0:
            limit_number_rows_select_stmt = ' limit {}'.format(str(number_rows_select))
    return limit_number_rows_select_stmt


def run_query(col_names_to_retrieve, cond_statement, limit_number_rows_select_stmt, cell_table, curs):
    curs.execute('select {} from {}{} {}'.format(col_names_to_retrieve, cell_table, cond_statement,
                                                 limit_number_rows_select_stmt))
    if curs is not None:
        return curs.fetchall()
    else:
        return []


def get_col_names_from_cells_assays(col_names, cells, assays, col_names_from_db):
    cells_to_report = []
    if len(cells) > 0:
        if cells[0] == 'all':
            for c in col_names_from_db:
                if '___' in c:
                    cells_to_report.append(c.split('___')[0])
        else:
            cells_to_report = cells
    for c in cells_to_report:
        for a in assays:
            if a == 'all':
                for ca in col_names_from_db:
                    if c in ca:  # if the cell name was in the col_names
                        if 'IMR__90'.lower() in c and 'IMR__90___RepliDomain'.lower() not in ca:  # for this cell just use the the replication timing (Hard coding :( )
                            continue
                        if ca not in col_names:
                            col_names.append(ca)
                break
            elif c + "___" + a in col_names_from_db:
                if c + "___" + a not in col_names:
                    col_names.append(c + "___" + a)
    return ','.join(col_names)


def get_cell_info_for_regions(regions_input_file, db_name='testregmotifs', cell_table='motifs',
                              cells=['HepG2'], assays=['all'],
                              col_names=['chr', 'motifstart', 'motifend', 'name', 'score', 'pval', 'strand'],
                              number_rows_select='all',
                              sep='\t', report_cols_from_file=True, cols_indices_to_report_from_file=[9],
                              cols_names_to_report_from_file=['Activity_Score'], df_output_file='',
                              region_name_index=6, region_strand_index=4, region_score_index=4, motif_score_index=4,
                              max_number_motifs_to_report=1,
                              min_dist_from_region_start=True, min_dist_from_region_center=True, max_motif_score=True,
                              max_region_score=True
                              ):
    conn = open_connection(db_name)
    curs = conn.cursor()
    cells = updateColNames(cells)
    assays = updateColNames(assays)
    col_names_from_db = get_col_names_from_table(cell_table, conn)
    limit_stmt = get_limit_smt(number_rows_select=number_rows_select)
    col_names_to_retrieve = get_col_names_from_cells_assays(col_names, cells, assays, col_names_from_db)
    results_regions = {}
    reported_col_names = col_names_to_retrieve.split(',')
    reported_col_names.extend(cols_names_to_report_from_file)
    if os.path.exists(df_output_file):
        close_connection(conn)
        return pd.read_pickle(df_output_file)

    number_lines_processed = 0
    with open(regions_input_file, 'r') as regions_infile:
        line = regions_infile.readline()
        print('Get motif info for', cells)
        while line:
            sline = line.strip().split(sep)
            conds = []
            motif_info_to_query = [sline[0], sline[1], sline[2]]  # get it from the file

            if len(motif_info_to_query) >= 3:
                conds.append('(chr={} and (posrange && int4range({},{})))'.format(
                    int(sline[0].replace('X', '23').replace('Y', '24').replace('M', '25').replace('chr', '')),
                    int(sline[1]), int(sline[2]) + 1))
            cond_stmt = get_cond_stmt(conds)
            rows = run_query(col_names_to_retrieve, cond_stmt, limit_stmt, cell_table, curs)
            cols_to_report_from_file = []
            for c in cols_indices_to_report_from_file:
                cols_to_report_from_file.append(sline[c])
            for row in rows:
                lrow = list(row)
                if 'NaN' in lrow or 'nan' in lrow or (str(float('nan')) in [str(x) for x in lrow]):
                    continue
                lrow.extend(cols_to_report_from_file)
                if sline[region_name_index] not in results_regions.keys():
                    results_regions[sline[region_name_index]] = [lrow]
                elif max_number_motifs_to_report == 'all':
                    results_regions[sline[region_name_index]].append(lrow)
                elif len(results_regions[sline[region_name_index]]) < max_number_motifs_to_report:
                    results_regions[sline[region_name_index]].append(lrow)
                else:  # if more than the require motifs was already reported then try to replace reported ones with better (how?) motifs.
                    for i in range(0, len(results_regions[sline[region_name_index]])):
                        if min_dist_from_region_start:
                            if sline[region_strand_index] == '-':
                                if abs(int(lrow[1]) - int(sline[1])) < abs(
                                        int(results_regions[sline[region_name_index]][i][1]) - int(sline[1])):
                                    results_regions[sline[region_name_index]][i] = lrow
                                    break
                            else:
                                if abs(int(lrow[2]) - int(sline[2])) < abs(
                                        int(results_regions[sline[region_name_index]][i][2]) - int(sline[2])):
                                    results_regions[sline[region_name_index]][i] = lrow
                                    break
                        elif min_dist_from_region_center:
                            if abs(int(lrow[1]) - ((int(sline[1]) + int(sline[1])) / 2)) < abs(
                                    int(results_regions[sline[region_name_index]][i][1]) - (
                                            (int(sline[1]) + int(sline[1])) / 2)):
                                results_regions[sline[region_name_index]][i] = lrow
                                break
                        elif max_motif_score:
                            if float(lrow[motif_score_index]) > results_regions[sline[region_name_index]][i][
                                motif_score_index]:
                                results_regions[sline[region_name_index]][i] = lrow
                                break
                        elif max_region_score:
                            if float(sline[region_score_index]) > results_regions[sline[region_name_index]][i][
                                region_score_index + 6]:
                                results_regions[sline[region_name_index]][i] = lrow
                                break
            line = regions_infile.readline()
            number_lines_processed += 1
            if number_lines_processed % 500 == 0:
                print('number_lines_processed: ', number_lines_processed)

    close_connection(conn)
    results_lines = []
    for region in results_regions.keys():
        for motif_region in results_regions[region]:
            results_lines.append(motif_region)
    results_df = pd.DataFrame(data=results_lines, columns=reported_col_names)
    results_df.to_pickle(df_output_file)

    return results_df


def get_promoters_of_unactive_genes(genes_input_file, proms_of_unactive_genes_output_file, unactive_expr_thresh=0,
                                    gene_expression_index=-1, strand_index=5, sep='\t', num_bp_for_prom=1000):
    if os.path.exists(proms_of_unactive_genes_output_file):
        return proms_of_unactive_genes_output_file
    with open(genes_input_file, 'r') as genes_infile, open(proms_of_unactive_genes_output_file, 'w') as outfile:
        print('Getting promoters (size={}) of KNOWN protein_coding genes that have expr <= {} in {}'.format(
            num_bp_for_prom,
            unactive_expr_thresh,
            genes_input_file))
        l = genes_infile.readline()
        while l:
            sl = l.strip().split(sep)
            if float(sl[gene_expression_index]) <= unactive_expr_thresh and sl[7] == 'KNOWN' and \
                    sl[8] == 'protein_coding':  # 'pseudogene' not in sl[8]:
                if sl[strand_index] == "-":
                    outfile.write(
                        sl[0] + sep + sl[2] + sep + str(int(sl[2]) + num_bp_for_prom) + sep + sep.join(sl[3::]) + '\n')
                else:
                    outfile.write(
                        sl[0] + sep + str(int(sl[1]) - num_bp_for_prom) + sep + sl[1] + sep + sep.join(sl[3::]) + '\n')
            l = genes_infile.readline()
    return proms_of_unactive_genes_output_file


def get_unactive_motifs(sys_args, db_name, prom_unactive_genes_start_index_params, col_names):
    genes_input_file = sys_args[prom_unactive_genes_start_index_params]
    proms_of_unactive_genes_output_file = sys_args[prom_unactive_genes_start_index_params + 1]
    prom_df_output_file = sys_args[prom_unactive_genes_start_index_params + 2]
    cells_to_extract_info_from_for_prom_unactive = sys_args[prom_unactive_genes_start_index_params + 3].split(',')
    proms_unactive_genes = get_promoters_of_unactive_genes(genes_input_file=genes_input_file,
                                                           proms_of_unactive_genes_output_file=proms_of_unactive_genes_output_file,
                                                           unactive_expr_thresh=0, gene_expression_index=-1,
                                                           strand_index=5, sep='\t', num_bp_for_prom=1000)
    prom_results_df = get_cell_info_for_regions(proms_unactive_genes, db_name=db_name,
                                                cells=cells_to_extract_info_from_for_prom_unactive, assays=['all'],
                                                df_output_file=prom_df_output_file, col_names=col_names)
    prom_results_df.to_csv(prom_df_output_file + '.tsv', sep='\t')


if __name__ == '__main__':
    db_name = "funmotifsdb"
    args = ['../datafiles/GeneExp/ENCODE/HepG2.bed', '../datafiles/TrainingSets/HepG2_unactive_proms.bed',
            '../datafiles/TrainingSets/HepG2_unactive_proms_motifs.df', 'HepG2,Liver']
    prom_unactive_genes_start_index_params = 0
    col_names = ['chr', 'motifstart', 'motifend', 'name', 'score', 'pval', 'strand']
    get_unactive_motifs(args, db_name, prom_unactive_genes_start_index_params, col_names)
    print("Done.")
    # TODO: write unit test