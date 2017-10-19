'''
Created on 15 Oct 2017

@author: husensofteng
'''
from multiprocessing import Pool
import DBUtilities

def generate_PFM_table(PFMs_dict, PFM_table_name,
                       db_name, db_user_name, db_host_name,
                       cols,
                       cols_names):
    '''
    Store the PFM of each motif in PFM_table
    each contains: motif_name, position in the motif, allele of the position, frequency of this allele in this position
    '''
    rows_to_insert = []
    for motif in PFMs_dict.keys():
        for pos in range(0, len(PFMs_dict[motif])):
            PFMs_dict[motif][0]
            for allele in PFMs_dict[motif][pos].keys():
                rows_to_insert.append([motif, pos+1, allele, PFMs_dict[motif][pos][allele]])
    
    value_marks = ['%s' for i in range(0, len(rows_to_insert[0]))]
    conn = DBUtilities.open_connection(db_name, db_user_name, db_host_name)
    curs = conn.cursor()
    curs.execute('drop table if exists {0};'.format(PFM_table_name))
    curs.execute('create table if not exists {0} ({1});'.format(PFM_table_name, ','.join(cols)))
    conn.commit()
    curs.executemany('insert into {} ({}) values({})'.format(PFM_table_name, ', '.join(cols_names), ', '.join(value_marks)), rows_to_insert)
    curs.execute('create index if not exists {0}name on {1} using btree(name);'.format(PFM_table_name, PFM_table_name))
    conn.commit()
    curs.close()
    conn.close()
    
 
def create_motifs_table(db_name, db_user_name, db_host_name, motifs_table, motif_cols, new_table_name):
    conn = DBUtilities.open_connection(db_name, db_user_name, db_host_name)
    curs = conn.cursor()
    
    curs.execute('drop table if exists {0};'.format(new_table_name))
    print 'create table if not exists {0} as (select {1} from {2});'.format(new_table_name, ','.join(motif_cols), motifs_table)
    curs.execute('create table if not exists {0} as (select {1} from {2});'.format(new_table_name, ','.join(motif_cols), motifs_table))
    curs.execute('create index if not exists {0}mid on {1} using btree(mid);'.format(new_table_name, new_table_name))
    curs.execute('create index if not exists {0}posrange on {1} using gist(posrange);'.format(new_table_name, new_table_name))
    curs.execute('create index if not exists {0}tfname on {1} using btree(name);'.format(new_table_name, new_table_name))
    curs.execute('create index if not exists {0}chr on {1} using btree(chr);'.format(new_table_name, new_table_name))
    conn.commit()
    curs.close()
    conn.close()

def split_motifs_parallel(db_name, db_user_name, db_host_name, motifs_table, chr, motif_cols):
    conn = DBUtilities.open_connection(db_name, db_user_name, db_host_name)
    curs = conn.cursor()
    new_table_name = "chr"+str(chr)+"motifs"
    curs.execute('drop table if exists {0};'.format(new_table_name))
    print 'create table if not exists {0} as (select {1} from {2} where chr={3});'.format(new_table_name, ','.join(motif_cols), motifs_table, chr)
    curs.execute('create table if not exists {0} as (select {1} from {2} where chr={3});'.format(new_table_name, ','.join(motif_cols), motifs_table, chr))
    curs.execute('create index if not exists {0}mid on {1} using btree(mid);'.format("chr"+str(chr), new_table_name))
    curs.execute('create index if not exists {0}posrange on {1} using gist(posrange);'.format("chr"+str(chr), new_table_name))
    conn.commit()
    curs.close()
    conn.close()
    

def split_motifs_table_by_chr(db_name, db_user_name, db_host_name, 
                              motifs_table, 
                              motif_cols, 
                              chr_names):
    
    p = Pool()
    for chr in chr_names:
        p.apply_async(split_motifs_parallel, args = (db_name, db_user_name, db_host_name, motifs_table, chr, motif_cols))
        #split_motifs_parallel(db_name, db_user_name, db_host_name, motifs_table, chr, motif_cols)
    p.close()
    p.join()
    return

                       
def motif_names_table(db_name, db_user_name, db_host_name, 
                      motifs_table,
                      motif_names_table):
    conn = DBUtilities.open_connection(db_name, db_user_name, db_host_name)
    curs = conn.cursor()
    curs.execute('drop table if exists {0};'.format(motif_names_table))
    print 'create table if not exists {0} as (select distinct(name) from {1});'.format(motif_names_table, motifs_table)
    curs.execute('create table if not exists {0} as (select distinct(name) from {1});'.format(motif_names_table, motifs_table))
    conn.commit()
    curs.close()
    
