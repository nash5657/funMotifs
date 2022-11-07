'''
Created on Dec 11, 2016

@author: husensofteng
'''
import sys, os
import pandas as pd
from DBUtilities import open_connection, close_connection, get_col_names_from_table


def updateColNames(col_names=[]):
    for i in range(0,len(col_names)):
        col_names[i] = ('_'.join((col_names[i].replace('(','').replace(')','').replace('-','__')).split())).lower()
    return col_names

def get_cond_stmt(conds = []):
    if len(conds)>0:
        return ' where {} '.format(' and '.join(conds))
    else:
        return ''
    
def get_limit_smt(number_rows_select='all'):
    limit_number_rows_select_stmt = ""
    if number_rows_select!="all":
        if number_rows_select>0:
            limit_number_rows_select_stmt = ' limit {}'.format(str(number_rows_select))
    return limit_number_rows_select_stmt

def run_query(col_names_to_retrieve, cond_statement, limit_number_rows_select_stmt, cell_table, curs):
    
    curs.execute('select {} from {}{} {}'.format(col_names_to_retrieve, cell_table, cond_statement, limit_number_rows_select_stmt))
    if curs is not None:
        return curs.fetchall()
    else:
        return []

def get_col_names_from_cells_assays(col_names, cells, assays, col_names_from_db):
    cells_to_report = []
    if len(cells)>0:
        if cells[0]=='all':
            for c in col_names_from_db:
                if '___' in c:
                    cells_to_report.append(c.split('___')[0])
        else:
            cells_to_report = cells
    for c in cells_to_report:
        for a in assays:
            if a == 'all':
                for ca in col_names_from_db:
                    if c in ca:#if the cell name was in the col_names
                        if 'IMR__90'.lower() in c and 'IMR__90___RepliDomain'.lower() not in ca:#for this cell just use the the replication timing (Hard coding :( )
                            continue
                        if ca not in col_names:
                            col_names.append(ca)
                break
            elif c+"___"+a in col_names_from_db:
                if c+"___"+a not in col_names:
                    col_names.append(c+"___"+a)
    return ','.join(col_names)
    
def get_cell_info_for_motifs(motifs_input_file, db_name = 'funmotifsdb', cell_table='motifs', cells = ['HepG2'], assays = ['all'],
                             col_names = ['chr', 'motifstart', 'motifend', 'name', 'score', 'pval', 'strand'], motif_tf_names_to_consider=[], number_rows_select='all',
                             sep='\t', report_cols_from_file=True, cols_indices_to_report_from_file=[7], cols_names_to_report_from_file=['Activity_Score'], df_output_file=''):
    
    # TODO: change db_user_name
    conn = open_connection(db_name, db_user_name="markmzr")
    if os.path.exists(df_output_file):
        close_connection(conn)
        return pd.read_pickle(df_output_file)
    
    curs = conn.cursor()
    cells = updateColNames(cells)
    assays = updateColNames(assays)
    col_names_from_db = get_col_names_from_table(cell_table, conn)
    limit_stmt = get_limit_smt(number_rows_select=number_rows_select) 
    col_names_to_retrieve = get_col_names_from_cells_assays(col_names, cells, assays, col_names_from_db)
    results_regions = {}
    reported_col_names = col_names_to_retrieve.split(',')
    reported_col_names.extend(cols_names_to_report_from_file)
    
    number_lines_processed = 0
    with open(motifs_input_file, 'r') as motif_infile:
        line = motif_infile.readline()
        print(('Get motif info for', cells))
        while line:
            sline = line.strip().split(sep)
            sline[7] = float(sline[7])
            #only take active motifs
            if abs(sline[7])<0.5:
                line = motif_infile.readline()
                continue
            conds = []
            motif_info_to_query = [sline[0], sline[1], sline[2], sline[3]]#get it from the file
            
            if len(motif_info_to_query)>=4:
                conds.append('(' + "chr={} and (posrange = int4range({}, {})) and motifstart={} and motifend={} and name='{}'".format(int(motif_info_to_query[0].replace('X', '23').replace('Y', '24').replace('M', '25').replace('chr', '')), 
                                                                                                    int(motif_info_to_query[1]), int(motif_info_to_query[2])+1, 
                                                                                                    int(motif_info_to_query[1]), int(motif_info_to_query[2]), motif_info_to_query[3]) + ')')
            
            cond_stmt = get_cond_stmt(conds)
            rows = run_query(col_names_to_retrieve, cond_stmt, limit_stmt, cell_table, curs)
            cols_to_report_from_file = []
            for c in cols_indices_to_report_from_file:
                cols_to_report_from_file.append(sline[c])
            for row in rows:
                lrow = list(row)
                MPRA_score_index = len(lrow) + 1
                if 'NaN' in lrow or 'nan' in lrow or str(float('nan')) in [str(x) for x in lrow]:#postgres converts NaNs to float('nan') by default
                        continue
                lrow.extend(cols_to_report_from_file)
                if lrow[-1] not in list(results_regions.keys()): 
                    results_regions[sline[-1]] = [lrow]
                elif len(results_regions[sline[-1]])<1:
                    results_regions[sline[-1]].append(lrow)
                else:
                    for i in range(0, len(results_regions[sline[-1]])):
                        if abs(lrow[MPRA_score_index]) > abs(results_regions[sline[-1]][i][MPRA_score_index]):
                            results_regions[sline[-1]][i] = lrow
            line = motif_infile.readline()
            number_lines_processed+=1
            if number_lines_processed % 500 == 0:
                print(('number_lines_processed: ', number_lines_processed))
    close_connection(conn)
    results_lines = []
    for region in list(results_regions.keys()):
        for motif_region in results_regions[region]:
            results_lines.append(motif_region)
    results_df = pd.DataFrame(data=results_lines, columns=reported_col_names)
    results_df.to_pickle(df_output_file)
    
    return results_df


def get_cell_info_for_regions(regions_input_file, db_name = 'testregmotifs', cell_table='cell_table',
                              cells = ['HepG2'], assays = ['all'], 
                              col_names = ['chr', 'motifstart', 'motifend', 'name', 'score', 'pval', 'strand'], number_rows_select='all',
                              sep='\t', report_cols_from_file=True, cols_indices_to_report_from_file=[9], cols_names_to_report_from_file=['Activity_Score'], df_output_file='',
                              region_name_index = 6, region_strand_index = 4, region_score_index = 4, motif_score_index = 4, max_number_motifs_to_report = 1,
                              min_dist_from_region_start = True, min_dist_from_region_center = True, max_motif_score = True, max_region_score = True
                              ):
    # TODO: db_user_name as argument

    conn = open_connection(db_name, db_user_name='markmzr')
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
        print(('Get motif info for', cells))
        while line:
            sline = line.strip().split(sep)
            conds = []
            motif_info_to_query = [sline[0], sline[1], sline[2]]#get it from the file
            
            if len(motif_info_to_query)>=3:
                conds.append('(chr={} and (posrange && int4range({},{})))'.format(int(sline[0].replace('X', '23').replace('Y', '24').replace('M', '25').replace('chr', '')), 
                                                                                  int(sline[1]), int(sline[2])+1))
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
                if sline[region_name_index] not in list(results_regions.keys()): 
                    results_regions[sline[region_name_index]] = [lrow]
                elif max_number_motifs_to_report=='all':
                    results_regions[sline[region_name_index]].append(lrow)
                elif len(results_regions[sline[region_name_index]])<max_number_motifs_to_report:
                    results_regions[sline[region_name_index]].append(lrow)
                else:#if more than the require motifs was already reported then try to replace reported ones with better (how?) motifs.
                    for i in range(0, len(results_regions[sline[region_name_index]])):
                        if min_dist_from_region_start:
                            if sline[region_strand_index]=='-':
                                if abs(int(lrow[1])-int(sline[1])) < abs(int(results_regions[sline[region_name_index]][i][1]) - int(sline[1])):
                                    results_regions[sline[region_name_index]][i] = lrow
                                    break
                            else:
                                if abs(int(lrow[2])-int(sline[2])) < abs(int(results_regions[sline[region_name_index]][i][2]) - int(sline[2])):
                                    results_regions[sline[region_name_index]][i] = lrow
                                    break
                        elif min_dist_from_region_center:
                            if abs(int(lrow[1])-((int(sline[1])+int(sline[1]))/2)) < abs(int(results_regions[sline[region_name_index]][i][1]) - ((int(sline[1])+int(sline[1]))/2)):
                                results_regions[sline[region_name_index]][i] = lrow
                                break
                        elif max_motif_score:
                            if float(lrow[motif_score_index]) > results_regions[sline[region_name_index]][i][motif_score_index]:
                                results_regions[sline[region_name_index]][i] = lrow
                                break
                        elif max_region_score:
                            if float(sline[region_score_index]) > results_regions[sline[region_name_index]][i][region_score_index+6]:
                                results_regions[sline[region_name_index]][i] = lrow
                                break
            line = regions_infile.readline()
            number_lines_processed+=1
            if number_lines_processed % 500 == 0:
                print(('number_lines_processed: ', number_lines_processed))
            
    close_connection(conn)
    results_lines = []
    for region in list(results_regions.keys()):
        for motif_region in results_regions[region]:
            results_lines.append(motif_region)
    results_df = pd.DataFrame(data=results_lines, columns=reported_col_names)
    results_df.to_pickle(df_output_file)
    
    return results_df
    
