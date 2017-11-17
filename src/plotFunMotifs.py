'''
Created on 17 Nov 2017

@author: husensofteng
'''
import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import tight_layout
import matplotlib.pyplot as plt

import sys
import pandas as pd
import seaborn as sns
import psycopg2
sns.set_style("white")
from pybedtools import BedTool
#plt.style.use('ggplot')
#sns.set_context("paper")#talk
#plt.style.use('seaborn-ticks')

params = {'-sep': '\t', 
          '-cols_to_retrieve':'chr, motifstart, motifend, strand, name, score, pval, fscore, chromhmm, contactingdomain, dnase__seq, fantom, loopdomain, numothertfbinding, othertfbinding, replidomain, tfbinding, tfexpr', '-number_rows_select':'all',
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
    
def plot_motif_freq(tf_name, tissue_table, motifs_table, min_fscore):
    
    conn = open_connection()
    curs = conn.cursor()
    
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
    curs = conn.cursor()
    
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
    curs = conn.cursor()
    
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
    curs = conn.cursor()
    
    stmt_all = "select upper(split_part(name,'_', 1)), count(name) as freq from {motifs},{tissue} where {motifs}.mid={tissue}.mid and ({tissue}.fscore>{min_fscore} or (tfbinding>0 and tfbinding!='NaN')) group by name order by freq desc".format(
        motifs=motifs_table, tissue=tissue_table, min_fscore=min_fscore)
    print stmt_all
    curs.execute(stmt_all)
    scores_all = curs.fetchall()
    curs.close()
    df = pd.DataFrame(scores_all, columns=['TFs', 'Frequency'])
    df['Tissue']=[tissue_table for i in range(0,len(df))]
    return df

def get_funmotifs(tissue_tables):
    cols = ['chr', 'motifstart', 'motifend', 'name', 'fscore', 'chromhmm', 'dnase__seq', 'contactingdomain', 'replidomain','tfbinding']
    
    conn = open_connection()
    curs = conn.cursor()
    motifs_table = 'motifs'
    for tissue_table in sorted(tissue_tables):
        print tissue_table
        query_stmt = "select {cols} from {motifs},{tissue} where {motifs}.mid={tissue}.mid and tfexpr>0 and ((fscore>2.0 and dnase__seq>0.0 and dnase__seq!='NaN' and tfbinding>0) or (tfbinding>0 and tfbinding!='NaN'))".format(
            cols=','.join(cols), motifs=motifs_table, tissue=tissue_table)
        print query_stmt
        curs.execute(query_stmt)
        query_results = curs.fetchall()
        df = pd.DataFrame(query_results, columns=cols)
        df.to_csv(tissue_table, sep='\t', index=False)
        df_bed = BedTool.from_dataframe(df).sort().merge(c=[4,5,6,7,8],o=['distinct','max', 'distinct', 'max','max', 'distinct', 'max'])
        df = BedTool.to_dataframe(df_bed, names=cols)
        df.to_csv(tissue_table+'_merged', sep='\t', index=False)
        print df.head()
    curs.close()
    conn.close()
    
if __name__ == '__main__':
    
    if len(sys.argv)<=0:
        print "Usage: python plotFunMotifs -plot"
        sys.exit(0)
    get_params(sys.argv[1:], params_without_value=[])
    tissue_tables=['blood', 'brain', 'breast','cervix', 'colon', 'esophagus', 'kidney', 'liver', 'lung', 'myeloid', 'pancreas', 'prostate', 'skin', 'stomach', 'uterus']
        
    get_funmotifs(tissue_tables)
    
    if '-plot' in params.keys():
        motifs_table='motifs'
        min_fscore = 2.5
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
            