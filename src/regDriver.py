'''
Created on 21 Oct 2017

@author: husensofteng
'''
import sys, os
import pandas as pd
import numpy as np
import string
import psycopg2
from multiprocessing import Pool
from psycopg2.extras import DictCursor
import time

params = {'-sep': '\t', '-cols_to_retrieve':'fscore', '-number_rows_select':'all',
          '-restart_conn_after_n_queries':100000, '-variant':True, '-region':True,
          '-chr':0, '-start':1, '-end':2, '-ref':3, '-alt':4, 
          '-db_name':'regmotifsdbtest', '-db_host':'localhost', '-db_port':5432, '-db_user':'huum', '-db_password':''}
    
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


def run_query(cols_to_retrieve, from_tabes, cond_statement, conn, n):
    curs = conn.cursor(name = "countcurs"+n, cursor_factory=DictCursor)
    stmt = 'select {} from {}{} {}'.format(cols_to_retrieve, from_tabes, cond_statement, get_limit_smt())
    curs.execute(stmt)
    if curs is not None:
        return curs.fetchall()
        curs.close()
    else:
        curs.close()
        return []

def read_infile():
    conn = open_connection()
    
    number_lines_processed = 0
    t = time.time()
    with open(params['-f'], 'r') as infile:
        line = infile.readline()
        while line:
            sline = line.strip().split(params['-sep'])
            if (line.startswith('#') or line.startswith('//') or len(sline)<3):
                line = infile.readline()
                continue
            if params['-variant']:#the input is variant
                try:
                    if ( #check if the number of ref/alt alleles match the variant length
                        (int(float(sline[params['-end']])) - int(float(sline[params['-start']])) + 1 != len(sline[params['-ref']]) and 
                         sline[params['-ref']]!='-' and sline[params['-alt']]!='-')):#skip mis appropriate lines
                            print 'Warning -- skipped line: the variant length does not match the ref/alt length', line
                            line = infile.readline()
                            continue
                except IndexError:
                    print 'Warning -- line is not a variant (fewer than 5 columns (chr,start,end,ref,alt) detected): ', line
                    params['-variant'] = False
                    
            updated_chr = sline[params['-chr']].replace('X', '23').replace('Y', '24').replace('MT','25').replace('M','25')
            chr_table = updated_chr+'motifs'
            if not updated_chr.startswith('chr'):
                chr_table = 'chr'+updated_chr+'motifs'
            cond_statement = (" where (posrange && int4range({start},{end},'[]')) and ({tissue_table}.mid={motif_table}.mid)".format(
                start=int(float(sline[params['-start']])), 
                end=int(float(sline[params['-end']])), 
                tissue_table=params['-tissue'], 
                motif_table=chr_table))
            
            rows = run_query(params['-cols_to_retrieve'], params['-tissue']+',' + chr_table, cond_statement, conn, str(number_lines_processed))
            print rows
            line = infile.readline()
            
            number_lines_processed+=1
            if number_lines_processed % int(params['-restart_conn_after_n_queries']) == 0:
                print '{} Lines are processed from {}'.format(number_lines_processed, params['-f'])
                print time.time()-t
                t = time.time()
                conn.close()
                conn = open_connection()
    return number_lines_processed
    
    
def get_motif_breaking_score(TF_motif_weights_dict, motif_name, motif_strand, motif_start, motif_end, mut_start, mut_end, ref_allele, alt_allele):
    
    if motif_strand=='-':
        ref_allele = ref_allele.translate(string.maketrans('ACGT','TGCA'))
        alt_allele = alt_allele.translate(string.maketrans('ACGT','TGCA'))
    
    breaking_score = 0.0
    breaking_score_cumulative = 0.0
    mut_sig = ""
    motif_mut_pos_start = 0
    motif_mut_pos_end = 0
    motif_length = motif_end-motif_start
    if mut_start >= motif_start and mut_end <=motif_end:#motif contains the mutation
        if motif_strand=='+':
            motif_mut_pos_start = mut_start-motif_start
            motif_mut_pos_end = mut_end-motif_start
        else:
            motif_mut_pos_start = motif_end-mut_end
            motif_mut_pos_end = motif_end-mut_start
    elif mut_start < motif_start and (mut_end >=motif_start and mut_end <=motif_end):#mut stretches to the left of the motif
        bp_to_strip = motif_start-mut_start
        if motif_strand == '+':
            motif_mut_pos_start = 0
            motif_mut_pos_end = mut_end-motif_start
        else:
            motif_mut_pos_start = motif_end-mut_end
            motif_mut_pos_end = motif_length
        
        if not ref_allele == '-':#if it is not insertion
            ref_allele = ref_allele[bp_to_strip:]
        if not alt_allele == '-' and not ref_allele == '-':#if it is not deletion nor insertion (don't touch insertions)
            alt_allele = alt_allele[bp_to_strip:]
            
            
    elif (mut_start >= motif_start and mut_start <= motif_end) and mut_end >motif_end:#mut stretches to the right of the motif
        if not ref_allele == '-':
            bp_to_strip = len(ref_allele)-(mut_end-motif_end)
            ref_allele = ref_allele[:bp_to_strip]
        if not alt_allele == '-' and not ref_allele == '-':
            #bp_to_strip = len(ref_allele)-(mut_end-motif_end)
            alt_allele = alt_allele[:bp_to_strip]
        
        if motif_strand=='+':
            motif_mut_pos_start = mut_start-motif_start
            motif_mut_pos_end = motif_length
        else:
            motif_mut_pos_start = 0
            motif_mut_pos_end = motif_end-mut_start
        
    elif mut_start < motif_start and mut_end > motif_end:#motif contains the mutation
        motif_mut_pos_start = 0
        motif_mut_pos_end = motif_length
        bp_to_strip = motif_start-mut_start
        if not ref_allele == '-':
            ref_allele = ref_allele[bp_to_strip:bp_to_strip+motif_length+1]
        if not alt_allele == '-' and not ref_allele=='-':
            alt_allele = alt_allele[bp_to_strip:bp_to_strip+motif_length+1]
    
    '''print TF_motif_weights_dict[motif_name]
    print motif_name, len(TF_motif_weights_dict[motif_name]), motif_strand
    print motif_start, motif_end
    print mut_start, mut_end
    print motif_mut_pos_start, motif_mut_pos_end
    print breaking_score
    print ref_allele, '>', alt_allele
    print mut_sig
    '''
    if ref_allele == '-' or alt_allele == '-':
        breaking_score = 1.0
        breaking_score_cumulative = (motif_mut_pos_end-motif_mut_pos_start)+1#number of deleted or inserted bps
        if breaking_score_cumulative>motif_length+1:#in cases where an insertion streches more than the motif then just count the number of bases in the motif and ignore the rest
            breaking_score_cumulative = motif_length+1
    else:
        for i, mut_pos in enumerate(range(motif_mut_pos_start, motif_mut_pos_end+1)):
            try:
                breaking_score += abs(TF_motif_weights_dict[motif_name][mut_pos][ref_allele[i]] - TF_motif_weights_dict[motif_name][mut_pos][alt_allele[i]])
                breaking_score_cumulative+=breaking_score
            except KeyError:
                continue
            #keep the breaking_score at max 1
        if breaking_score>=1.0:
            breaking_score = 1.0
    mut_sig = ref_allele+">"+alt_allele        
    
    
    motif_mut_pos = str(motif_mut_pos_start+1) + '-' + str(motif_mut_pos_end+1)
    if motif_mut_pos_start==motif_mut_pos_end:
        motif_mut_pos = str(motif_mut_pos_start+1)
    
    return breaking_score, breaking_score_cumulative, mut_sig, motif_mut_pos

if __name__ == '__main__':
    
    if len(sys.argv)<=0:
        print "Usage: python regDriver.py input_file [options]"
        sys.exit(0)
    get_params(sys.argv[1:], params_without_value=[])
    print params
    try:
        read_infile()
    except KeyError:
        print "No value was found for one or more of the arguments:\n", params
        print "Usage: python regDriver.py -f file_name -tissue tissue_name"
        sys.exit(0)
    
    