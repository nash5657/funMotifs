'''
Created on 21 Oct 2017

@author: husensofteng
'''
import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import tight_layout
import matplotlib.pyplot as plt

import sys, os
import pandas as pd
import numpy as np
import seaborn as sns
import string
import psycopg2
from multiprocessing import Pool
from psycopg2.extras import DictCursor
import time
#plt.style.use('ggplot')
sns.set_style("white")
#sns.set_context("paper")#talk
#plt.style.use('seaborn-ticks')

params = {'-sep': '\t', '-cols_to_retrieve':'chr, motifstart, motifend, strand, name, score, pval, fscore, chromhmm, contactingdomain, dnase__seq, fantom, loopdomain, numothertfbinding, othertfbinding, replidomain, tfbinding, tfexpr', '-number_rows_select':'all',
          '-restart_conn_after_n_queries':100000, '-variants':True, '-regions':True,
          '-chr':0, '-start':1, '-end':2, '-ref':3, '-alt':4, 
          '-db_name':'regmotifsdbtest', '-db_host':'localhost', '-db_port':5432, '-db_user':'huum', '-db_password':'',
          '-all_motifs':True, '-motifs_tfbining':False, '-max_score_motif':False, '-motifs_tfbinding_otherwise_max_score_motif':False,
          '-verbose': True}
    
def get_params(params_list, params_without_value):
    global params
    for i, arg in enumerate(params_list):#priority is for the command line
        if arg.startswith('-'): 
            if arg in params_without_value:
                params[arg] = True
            else:
                try:
                    v = params_list[i+1]
                    if v.lower()=='yes' or v.lower()=='true':
                        v=True
                    elif v.lower()=='no' or v.lower()=='false':
                        v=False
                    params[arg] =  v
                except IndexError:
                    print "no value is given for parameter: ", arg 
    return params

def open_connection():
    conn = psycopg2.connect("dbname={} user={} password={} host={} port={}".format(params['-db_name'], params['-db_user'], params['-db_password'], params['-db_host'], params['-db_port']))
    return conn

    
def get_col_names_from_table(table_name, conn):
    curs = conn.cursor()
    curs.execute("select * FROM {} limit 1".format(table_name))
    return [desc[0] for desc in curs.description]

def get_limit_smt():
    limit_number_rows_select_stmt = ""
    if params['-number_rows_select']!="all":
        if int(params['-number_rows_select'])>0:
            limit_number_rows_select_stmt = ' limit {}'.format(str(params['-number_rows_select']))
    return limit_number_rows_select_stmt


def run_query(cols_to_retrieve, from_tabes, cond_statement, order_by_stmt, conn, n):
    curs = conn.cursor(name = "countcurs"+n, cursor_factory=DictCursor)
    stmt = 'select {} from {}{}{} {}'.format(cols_to_retrieve, from_tabes, cond_statement, order_by_stmt, get_limit_smt())
    curs.execute(stmt)
    if curs is not None:
        return curs.fetchall()
        curs.close()
    else:
        curs.close()
        return []

def run_query_nocursorname(cols_to_retrieve, from_tabes, cond_statement, curs):
    #curs = conn.cursor()
    stmt = 'select {} from {}{} {}'.format(cols_to_retrieve, from_tabes, cond_statement, get_limit_smt())
    curs.execute(stmt)
    if curs is not None:
        return curs.fetchall()
        #curs.close()
    else:
        #curs.close()
        return []
    
def read_infile():
    conn = open_connection()
    curs_for_pfms = conn.cursor()
    number_lines_processed = 0
    t = time.time()
    with open(params['-f'], 'r') as infile, open(params['-f']+'_annotated.tsv', 'w') as outfile:
        line = infile.readline()
        cols_from_file = ['cols'+str(i) for i in range(0,len(line.strip().split(params['-sep'])))]
        cols_from_file.extend((params['-cols_to_retrieve'] + ',mutposition,entropy').split(','))
        outfile.write(params['-sep'].join(cols_from_file) + '\n')
        
        while line:
            sline = line.strip().split(params['-sep'])
            if (line.startswith('#') or line.startswith('//') or len(sline)<3):
                line = infile.readline()
                continue
            if params['-variants']:#the input is variant
                try:
                    if ( #check if the number of ref/alt alleles match the variant length
                        (int(float(sline[params['-end']])) - int(float(sline[params['-start']])) + 1 != len(sline[params['-ref']]) and 
                         sline[params['-ref']]!='-' and sline[params['-alt']]!='-' and
                         sline[params['-ref']]!='deletion' and sline[params['-alt']]!='insertion' and
                         sline[params['-ref']]!='del' and sline[params['-alt']]!='ins'
                         )):#skip mis appropriate lines
                            if params['-verbose']:
                                print 'Warning -- skipped line: the variant length does not match the ref/alt length', line
                            line = infile.readline()
                            continue
                except IndexError:
                    if params['-verbose']:
                        print 'Warning -- line is not a variant (fewer than 5 columns (chr,start,end,ref,alt) detected): ', line
                    params['-variants'] = False
                    
            updated_chr = sline[params['-chr']].replace('X', '23').replace('Y', '24').replace('MT','25').replace('M','25')
            chr_table = updated_chr+'motifs'
            if not updated_chr.startswith('chr'):
                chr_table = 'chr'+updated_chr+'motifs'
            cond_statement = (" where (posrange && int4range({start},{end},'[]')) and ({tissue_table}.mid={motif_table}.mid and {tissue_table}.tfexpr>0.0)".format(
                start=int(float(sline[params['-start']])), 
                end=int(float(sline[params['-end']])), 
                tissue_table=params['-tissue'], 
                motif_table=chr_table))
            #if params['-variants'] then also retreive the affinity change directly from the query (need for if and else in postgres)
            mutation_position_stmt = ''
            order_by_stmt = ' order by fscore '
            if params['-variants']:
                mutation_position_stmt = ", (CASE WHEN (UPPER(posrange * int4range({start}, {end})) - LOWER(posrange * int4range({start}, {end}))>1) THEN 100 ELSE (CASE when STRAND='-' THEN (motifend-{start})+1 ELSE ({start}-motifstart)+1 END) END) as mutposition ".format(start=int(float(sline[params['-start']])), end=int(float(sline[params['-end']])))
            rows = run_query(params['-cols_to_retrieve']+mutation_position_stmt, params['-tissue']+',' + chr_table, cond_statement, order_by_stmt, conn, str(number_lines_processed))
            #for each row get the entropy
            motifs_with_tfbinding = []
            all_motifs = []
            for row in rows:
                entropy = 0.0
                if (row['mutposition']==100 or
                    sline[params['-ref']]=='-' or sline[params['-ref']]=='deletion' or sline[params['-ref']]=='del' or
                    sline[params['-alt']]=='-'  or sline[params['-alt']]=='insertion' or sline[params['-alt']]=='ins'
                             ):
                    entropy=1
                else:
                    rows_pfms = run_query_nocursorname(cols_to_retrieve="(select freq from motifs_pfm where position={mutposition} and name = '{motif_name}' and allele='{ref_allele}') - (select freq from motifs_pfm where position={mutposition} and name = '{motif_name}' and allele='{alt_allele}')".format(
                        mutposition=row['mutposition'], motif_name=row['name'], ref_allele=sline[params['-ref']], alt_allele=sline[params['-alt']]), 
                                           from_tabes='motifs_pfm', cond_statement=" where position={mutposition} and name='{motif_name}' and allele='{ref_allele}'".format(
                                               mutposition=row['mutposition'], motif_name=row['name'], ref_allele=sline[params['-ref']]), curs=curs_for_pfms)
                    try:
                        entropy = float(rows_pfms[0][0])
                    except TypeError:
                        entropy = 'NA'
                        if params['-verbose']:
                            print 'Warning: ref/alt allele are not correctly given in: ' +  line
                        pass
                if row['numothertfbinding']<=0.0:
                    row['othertfbinding'] = "None"
                
                lrow=list(row)
                lrow.append(entropy)
                if params['-all_motifs']:
                    outfile.write(line.strip() + params['-sep'] + params['-sep'].join(str(x) for x in lrow) + '\n')
                else:
                    if float(row['tfbinding'])>0.0:
                        motifs_with_tfbinding.append(lrow)
                all_motifs.append(lrow)
            if not params['-all_motifs']:
                if params['-motifs_tfbining']:#only report motifs that have tfbinding
                    for motif_tfbinding in motifs_with_tfbinding:
                        outfile.write(line.strip() + params['-sep'] + params['-sep'].join(str(x) for x in motif_tfbinding) + '\n')
                elif params['-max_score_motif']:#don't care about tfbinding just give the motif with the maximum score
                    try:
                        max_motif = all_motifs[0]
                        outfile.write(line.strip() + params['-sep'] + params['-sep'].join(str(x) for x in max_motif) + '\n')
                    except IndexError:
                        pass
                elif params['-motifs_tfbinding_otherwise_max_score_motif']:#if there is any motif with tfbinding return it otherwise return the motif that has the maximum score
                    if len(motifs_with_tfbinding)>0:
                        for motif_tfbinding in motifs_with_tfbinding:
                            outfile.write(line.strip() + params['-sep'] + params['-sep'].join(str(x) for x in motif_tfbinding) + '\n')
                    else:
                        try:
                            max_motif = all_motifs[0]
                            outfile.write(line.strip() + params['-sep'] + params['-sep'].join(str(x) for x in max_motif) + '\n')
                        except IndexError:
                            pass
            line = infile.readline()
            
            number_lines_processed+=1
            if number_lines_processed % int(params['-restart_conn_after_n_queries']) == 0:
                print '{} Lines are processed from {}'.format(number_lines_processed, params['-f'])
                print time.time()-t
                t = time.time()
                conn.close()
                curs_for_pfms.close()
                conn = open_connection()
                curs_for_pfms = conn.cursor()
    return number_lines_processed    

def plot_motif_freq(tf_name, tissue_table, motifs_table, min_fscore):
    
    conn = open_connection()
    curs = conn.cursor()#cursor_factory=DictCursor)
    
    stmt_all = "select count({tissue}.mid) from {motifs},{tissue} where {motifs}.mid={tissue}.mid and {motifs}.name like '%{tf_name}%'".format(motifs=motifs_table, tissue=tissue_table, tf_name=tf_name)
    print stmt_all
    stmt_tfbinding = "select count({tissue}.mid) from {motifs},{tissue} where {motifs}.mid={tissue}.mid and {motifs}.name like '%{tf_name}%' and ({tissue}.tfbinding>0 and {tissue}.tfbinding!='NaN')".format(motifs=motifs_table, tissue=tissue_table,tf_name=tf_name)
    stmt_dnase = "select count({tissue}.mid) from {motifs},{tissue} where {motifs}.mid={tissue}.mid and {motifs}.name like '%{tf_name}%' and ({tissue}.dnase__seq>0 and {tissue}.dnase__seq!='NaN')".format(motifs=motifs_table, tissue=tissue_table, tf_name=tf_name)
    stmt_active = "select count({tissue}.mid) from {motifs},{tissue} where {motifs}.mid={tissue}.mid and {motifs}.name like '%{tf_name}%' and ((fscore>{min_fscore} and dnase__seq>0 and dnase__seq!='NaN') or (tfbinding>0 and {tissue}.tfbinding!='NaN'))".format(motifs=motifs_table, tissue=tissue_table, tf_name=tf_name, min_fscore=min_fscore)
    curs.execute(stmt_all)
    motifs_all = curs.fetchall()
    curs.execute(stmt_tfbinding)
    tfbinding = curs.fetchall()
    curs.execute(stmt_dnase)
    dnase = curs.fetchall()
    curs.execute(stmt_active)
    active = curs.fetchall()
    curs.close()
    return [[tf_name, tissue_table, 'All Motifs', int(motifs_all[0][0])], 
            [tf_name, tissue_table, 'DHSs', int(dnase[0][0])], 
            [tf_name, tissue_table, 'Matching TFBSs', int(tfbinding[0][0])], 
            [tf_name, tissue_table, 'Active Motifs', int(active[0][0])]]

def plot_fscore(tf_name, tissue_table, motifs_table, tissue_names, fig_name):
    
    conn = open_connection()
    curs = conn.cursor()#cursor_factory=DictCursor)
    
    stmt_all = "select {tissue_names} from {motifs},{tissue} where {motifs}.mid={tissue}.mid and {motifs}.name like '%{tf_name}%'".format(
        tissue_names=','.join(sorted(tissue_names)), motifs=motifs_table, tissue=tissue_table, tf_name=tf_name)
    print stmt_all
    curs.execute(stmt_all)
    scores_all = curs.fetchall()
    curs.close()
    df = pd.DataFrame(scores_all, columns=tissue_names)
    s = sns.boxplot(data=df, color='grey')
    ss = s.get_figure()
    ss.savefig(fig_name+'.pdf', bbox_inches='tight')
    ss.savefig(fig_name+'.svg', bbox_inches='tight')
    return

def plot_heatmap(min_fscore, motifs_table,tissue_table, fig_name, threshold_to_include_tf):
    conn = open_connection()
    curs = conn.cursor()#cursor_factory=DictCursor)
    
    stmt_all = "select chromhmm, upper(split_part(name,'_', 1)), count(name) from {motifs},{tissue} where {motifs}.mid={tissue}.mid and ({tissue}.fscore>{min_fscore} or (tfbinding>0 and tfbinding!='NaN')) group by chromhmm,name order by chromhmm".format(
        motifs=motifs_table, tissue=tissue_table, min_fscore=min_fscore)
    print stmt_all
    curs.execute(stmt_all)
    scores_all = curs.fetchall()
    curs.close()
    df = pd.DataFrame(scores_all, columns=['Chromatin States', 'TFs', 'Frequency'])
    df_pivot = df.pivot('Chromatin States', 'TFs', 'Frequency')
    df_pivot_filtered = pd.DataFrame()
    for c in df_pivot.columns:
            if df_pivot[c].sum()>threshold_to_include_tf:
                df_pivot_filtered[c] = df_pivot[c] 
    print df_pivot_filtered.head()
    if len(df_pivot_filtered)>0:
        s = sns.heatmap(data=df_pivot_filtered, square=True, cbar=True)
        ss = s.get_figure()
        ss.savefig(fig_name+'.pdf', bbox_inches='tight')
        ss.savefig(fig_name+'.svg', bbox_inches='tight')
    
def plot_scatter_plot(min_fscore, motifs_table, tissue_table):
    conn = open_connection()
    curs = conn.cursor()#cursor_factory=DictCursor)
    
    stmt_all = "select upper(split_part(name,'_', 1)), count(name) as freq from {motifs},{tissue} where {motifs}.mid={tissue}.mid and ({tissue}.fscore>{min_fscore} or (tfbinding>0 and tfbinding!='NaN')) group by name order by freq desc".format(
        motifs=motifs_table, tissue=tissue_table, min_fscore=min_fscore)
    print stmt_all
    curs.execute(stmt_all)
    scores_all = curs.fetchall()
    curs.close()
    df = pd.DataFrame(scores_all, columns=['TFs', 'Frequency'])
    df['Tissue']=[tissue_table for i in range(0,len(df))]
    return df

if __name__ == '__main__':
    
    if len(sys.argv)<=0:
        print "Usage: python regDriver.py input_file [options]"
        sys.exit(0)
    get_params(sys.argv[1:], params_without_value=[])
    if '-f' in params.keys():
        try:
            read_infile()
        except KeyError:
            print "No value was found for one or more of the arguments:\n", params
            print "Usage: python regDriver.py -f file_name -tissue tissue_name"
    
    if '-plot' in params.keys():
        motifs_table='motifs'
        min_fscore = 2.5
        tissue_tables=['blood', 'brain', 'breast','cervix', 'colon', 'esophagus', 'kidney', 'liver', 'lung', 'myeloid', 'pancreas', 'prostate', 'skin', 'stomach', 'uterus']
        tfs = ['CTCF', 'CEBPB', 'FOXA1', 'KFL14', 'HNF4A', 'MAFK']
        threshold_to_include_tf_in_heatmap = 20000
        sns.despine(right=True, top=True, bottom=False, left=False)
        if '-fig1' in params.keys():
            print 'plotting figure 1'
            for tissue_table in ['blood', 'liver']:
                fig = plt.figure()
                tfs_freq = []
                for tf in sorted(tfs):
                    tfs_freq.extend(plot_motif_freq(tf_name=tf, tissue_table = tissue_table, motifs_table = motifs_table, min_fscore = min_fscore))
                df = pd.DataFrame(tfs_freq, columns = ['TFs', 'Tissue', 'Activity', 'Frequency'])
                fig = plt.figure()
                s = sns.barplot(x='TFs', y='Frequency', hue='Activity', data=df, estimator=sum)
                ss = s.get_figure()
                ss.savefig('fig1_'+tissue_table+'.pdf', bbox_inches='tight')
                ss.savefig('fig1_'+tissue_table+'.svg', bbox_inches='tight')
               
        if '-fig2' in params.keys():
            print 'plotting figure 2'
            for tf in sorted(tfs):
                fig = plt.figure(figsize=(12,6))
                plot_fscore(tf_name=tf, tissue_table='all_tissues', motifs_table=motifs_table, tissue_names=tissue_tables, fig_name='fig2_'+tf)
        
        if '-fig3' in params.keys():
            print 'plotting figure 3'
            for tissue_table in ['liver', 'blood', 'breast']:
                fig = plt.figure(figsize=(12,6))
                plot_heatmap(min_fscore = min_fscore, motifs_table=motifs_table,tissue_table=tissue_table, fig_name='fig3_'+tissue_table, threshold_to_include_tf=threshold_to_include_tf_in_heatmap)
        
        if '-fig4' in params.keys():
            print 'plotting figure 4' 
            dfs = []
            for tissue_table in tissue_tables:
                dfs.append(plot_scatter_plot(min_fscore, motifs_table, tissue_table))
            
            all_dfs = pd.concat(dfs)
            fig = plt.figure(figsize=(12,8))
            s = sns.stripplot(x='Tissue', y='Frequency', data=all_dfs, jitter=True)
            s.set_ylim([0,80000])
            ss = s.get_figure()
            ss.savefig('fig4.pdf', bbox_inches='tight')
            ss.savefig('fig4.svg', bbox_inches='tight')