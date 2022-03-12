'''
Created on 15 Oct 2017

@author: husensofteng
'''
from multiprocessing import Pool
import DBUtilities

def insert_from_file(motif_cols_names, i_file, n, db_name, db_user_name, 
                     db_host_name, motif_table, 
                     thread_num=0):
    
    n_processed = 0
    conn = DBUtilities.open_connection(db_name, db_user_name, db_host_name)
    curs = conn.cursor()
    with open(i_file, 'r') as ifile:
        print(i_file)
        while True:
            #next_n_lines = list(islice(ifile,n))
            lines_as_lists = [l.strip().split('\t') for l in list(islice(ifile,n))]
            if not lines_as_lists:
                break
            n_processed+=len(lines_as_lists)#next_n_lines)
            #lines_as_lists = [l.strip().split('\t') for l in next_n_lines]
            try:
                value_marks = ['%s' for i in range(0, len(lines_as_lists[0]))]
                try:
                    curs.executemany('insert into {} ({}) values({})'.format(motif_table, ', '.join(motif_cols_names), ', '.join(value_marks)), lines_as_lists)
                    conn.commit()
                    print("Thread {} for ({}) has processed: {}".format(thread_num, i_file, n_processed))
                except:
                    print("Thread {} coundn't insert to DB ({})".format(thread_num, i_file))
            except (ValueError, IndexError):
                print("Empty file" + i_file + "in thread" + str(thread_num))
    curs.close()
    conn.close()
    return

def insert_into_db(motif_cols_names, db_name, db_user_name, db_host_name,  
                       motif_table, 
                       scored_motifs_overlapping_tracks_files_in, 
                       run_in_parallel_param,
                       number_processes_to_run_in_parallel
                       ):#, dir_to_import, keyword_to_check, header):
    log_file = open("log_file.txt", 'w')
    if run_in_parallel_param and len(scored_motifs_overlapping_tracks_files_in)>1:
        p = Pool(number_processes_to_run_in_parallel)
    thread_num=0
    for i_file in scored_motifs_overlapping_tracks_files_in:#[f for f in glob.glob('{}/*{}*'.format(dir_to_import, keyword_to_check))]
        if run_in_parallel_param and len(scored_motifs_overlapping_tracks_files_in)>1:
            thread_num+=1
            p.apply_async(insert_from_file, args=[motif_cols_names, i_file, 100000, db_name, db_user_name, db_host_name, motif_table,  thread_num])
        else:
            insert_from_file(motif_cols_names, i_file, 100000, db_name, db_user_name, db_host_name, motif_table)
        log_file.write(i_file+'\n')
    if run_in_parallel_param and len(scored_motifs_overlapping_tracks_files_in)>1:
        p.close()
        p.join()
    
    print("Data insertion into {} is done".format(cell_table))
    return


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
    return
    
 
def create_motifs_table(db_name, db_user_name, db_host_name, motifs_table, motif_cols, new_table_name):
    conn = DBUtilities.open_connection(db_name, db_user_name, db_host_name)
    curs = conn.cursor()
    
    curs.execute('drop table if exists {0};'.format(new_table_name))
    print('create table if not exists {0} as (select {1} from {2});'.format(new_table_name, ','.join(motif_cols), motifs_table))
    curs.execute('create table if not exists {0} as (select {1} from {2});'.format(new_table_name, ','.join(motif_cols), motifs_table))
    curs.execute('create index if not exists {0}mid on {1} using btree(mid);'.format(new_table_name, new_table_name))
    curs.execute('create index if not exists {0}posrange on {1} using gist(posrange);'.format(new_table_name, new_table_name))
    curs.execute('create index if not exists {0}tfname on {1} using btree(name);'.format(new_table_name, new_table_name))
    curs.execute('create index if not exists {0}chr on {1} using btree(chr);'.format(new_table_name, new_table_name))
    conn.commit()
    curs.close()
    conn.close()
    return


def create_motifs_table_from_file(db_name, db_user_name, db_host_name, scored_motifs_overlapping_tracks_files, motif_cols, motif_cols_names, new_table_name, run_in_parallel_param,
                       number_processes_to_run_in_parallel):
    
    #create database
    conn = DBUtilities.open_connection(db_name, db_user_name, db_host_name)
    curs = conn.cursor()
    
    curs.execute('drop table if exists {0};'.format(new_table_name))
    curs.execute('create table if not exists {0} ({1});'.format(new_table_name, ','.join(motif_cols)))
    
    conn.commit()
    curs.close()
    DBUtilities.close_connection(conn)
    
    print("Inserting data into: ", new_table_name)
    #retrive motif info from files
    for file_in in scored_motifs_overlapping_tracks_files:
        comm_cut_files ="cut -f1-{} {}| tail -n +2 > {}".format(len(motif_cols_names), file_in,  file_in+"_motifs")
        print(comm_divide_files)
        os.system(comm_divide_files)
        
        insert_into_db(motif_cols_names,
                       db_name = db_name, 
                       db_user_name=db_user_name, 
                       db_host_name=db_host_name, 
                       motif_table=new_table_name, 
                       scored_motifs_overlapping_tracks_files_in=file_in, 
                       run_in_parallel_param=run_in_parallel_param,
                       number_processes_to_run_in_parallel=number_processes_to_run_in_parallel)
    
        os.remove(file_in)
    
    
    
    
    conn = DBUtilities.open_connection(db_name, db_user_name, db_host_name)
    curs = conn.cursor()
    curs.execute('create index if not exists {0}mid on {1} using btree(mid);'.format(new_table_name, new_table_name))
    curs.execute('create index if not exists {0}posrange on {1} using gist(posrange);'.format(new_table_name, new_table_name))
    curs.execute('create index if not exists {0}tfname on {1} using btree(name);'.format(new_table_name, new_table_name))
    curs.execute('create index if not exists {0}chr on {1} using btree(chr);'.format(new_table_name, new_table_name))
    conn.commit()
    curs.close()
    conn.close()
    return
    
    
def split_motifs_parallel(db_name, db_user_name, db_host_name, motifs_table, chr, motif_cols):
    conn = DBUtilities.open_connection(db_name, db_user_name, db_host_name)
    curs = conn.cursor()
    new_table_name = "chr"+str(chr)+"motifs"
    curs.execute('drop table if exists {0};'.format(new_table_name))
    print('create table if not exists {0} as (select {1} from {2} where chr={3});'.format(new_table_name, ','.join(motif_cols), motifs_table, chr))
    curs.execute('create table if not exists {0} as (select {1} from {2} where chr={3});'.format(new_table_name, ','.join(motif_cols), motifs_table, chr))
    curs.execute('create index if not exists {0}mid on {1} using btree(mid);'.format("chr"+str(chr), new_table_name))
    curs.execute('create index if not exists {0}posrange on {1} using gist(posrange);'.format("chr"+str(chr), new_table_name))
    conn.commit()
    curs.close()
    conn.close()
    return
    

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
    print('create table if not exists {0} as (select distinct(name) from {1});'.format(motif_names_table, motifs_table))
    curs.execute('create table if not exists {0} as (select distinct(name) from {1});'.format(motif_names_table, motifs_table))
    conn.commit()
    curs.close()
    return
