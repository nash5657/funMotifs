'''
Created on 28 Sep 2017

@author: husensofteng
'''
import os,sys
from collections import Counter
from itertools import islice
from multiprocessing import Pool
import time
import math
import psycopg2
from psycopg2.extensions import ISOLATION_LEVEL_AUTOCOMMIT
from psycopg2.extras import DictCursor

def db_setup(cells_assays_dict, 
             assay_cells_datatypes, 
             db_name, cell_table, 
             db_user_name, db_host_name,
             motif_col_names, 
             motif_cols
             ):
    conn = ""
    curs = ""
    try:
        conn = psycopg2.connect("dbname={} user={} host={}".format(db_name, db_user_name, db_host_name))
        curs = conn.cursor()
        curs.execute("SELECT exists(SELECT 1 from pg_catalog.pg_database where datname = %s)", (db_name,))    
        if curs.fetchone()[0]:
            print "Successfully connected to DB: ", db_name
    
    except psycopg2.DatabaseError, e:
        print "Error %s" %e
        print "Creating DB: ", db_name
        con_postgres = psycopg2.connect(dbname='postgres', user=db_user_name, host=db_host_name)
        con_postgres.set_isolation_level(ISOLATION_LEVEL_AUTOCOMMIT)
        curs_postgres = con_postgres.cursor()
        curs_postgres.execute('CREATE DATABASE ' + db_name)
        con_postgres.commit()
        curs_postgres.close()
        con_postgres.close()
        
        conn = psycopg2.connect("dbname={} user={} host={}".format(db_name, db_user_name, db_host_name))
        curs = conn.cursor()
        curs.execute("SELECT exists(SELECT 1 from pg_catalog.pg_database where datname = %s)", (db_name,))    
        if curs.fetchone()[0]:
            print "Successfully created and connected to DB: ", db_name
    
    field_names = []
    col_names = []
    field_names.extend(motif_cols)
    col_names.extend(motif_col_names )
    for cell in sorted(cells_assays_dict.keys()):
        for assay in sorted(cells_assays_dict[cell].keys()):
            data_type = 'text'
            try:
                data_type = assay_cells_datatypes[assay]
            except KeyError:
                pass
            field_names.append('_'.join(((cell + "___" + assay).replace('(','').replace(')','')
                                         .replace('-','__')).split()) + " " + data_type)
            col_names.append('_'.join(((cell + "___" + assay).replace('(','').replace(')','')
                                         .replace('-','__')).split()))
    #curs.execute("DROP TABLE IF EXISTS {}".format(cell_table))
    
    create_table_stmt = "CREATE TABLE IF NOT EXISTS {} ({});".format(cell_table, ' ,'.join(field_names))
    curs.execute(create_table_stmt)
    conn.commit()
    conn.close()
    return col_names

def open_connection(db_name, db_user_name, db_host_name):
    conn = psycopg2.connect("dbname={} user={} host={}".format(db_name, db_user_name, db_host_name))
    #conn.row_factory = psycopg2.Row
    return conn

def close_connection(conn):
    conn.close()

def get_col_names_from_table(table_name, conn):
    curs = conn.cursor()
    curs.execute("select * FROM {} limit 1".format(table_name))
    return [desc[0] for desc in curs.description]
    
def insert_from_file(field_names, i_file, n, db_name, db_user_name, 
                     db_host_name, cell_table, header,
                     thread_num=0):
    
    n_processed = 0
    conn = open_connection(db_name, db_user_name, db_host_name)
    curs = conn.cursor()
    with open(i_file, 'r') as ifile:
        print i_file
        if header:#read an extra line to skip the first one
            list(islice(ifile, 1))
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
                    curs.executemany('insert into {} ({}) values({})'.format(cell_table, ', '.join(field_names), ', '.join(value_marks)), lines_as_lists)
                    conn.commit()
                    print "Thread {} for ({}) has processed: {}".format(thread_num, i_file, n_processed)
                except:
                    print "Thread {} coundn't insert to DB ({})".format(thread_num, i_file)
            except (ValueError, IndexError):
                print "Empty file" + i_file + "in thread" + str(thread_num)
    curs.close()
    conn.close()
    return

def insert_into_db(field_names, db_name, db_user_name, db_host_name,  
                       cell_table, 
                       scored_motifs_overlapping_tracks_files, 
                       header,
                       run_in_parallel_param,
                       number_processes_to_run_in_parallel
                       ):#, dir_to_import, keyword_to_check, header):
    
    if run_in_parallel_param and len(scored_motifs_overlapping_tracks_files)>1:
        p = Pool(number_processes_to_run_in_parallel)
    thread_num=0
    for i_file in scored_motifs_overlapping_tracks_files:#[f for f in glob.glob('{}/*{}*'.format(dir_to_import, keyword_to_check))]
        if run_in_parallel_param and len(scored_motifs_overlapping_tracks_files)>1:
            thread_num+=1
            p.apply_async(insert_from_file, args=[field_names, i_file, 1000000, db_name, db_user_name, db_host_name, cell_table, header, thread_num])
        else:
            insert_from_file(field_names, i_file, 1000000, db_name, db_user_name, db_host_name, cell_table, header)
    
    if run_in_parallel_param and len(scored_motifs_overlapping_tracks_files)>1:
        p.close()
        p.join()
    print "Data insertion into {} is done".format(cell_table)
    return

def get_tissue_cell_mappings(cell_assays, assay_names, 
                             tissue_cell_mappings_file, 
                             key_value_sep='=', 
                             values_sep=',', 
                             cell_assay_sepe=':',
                             motif_cols = ['posrange', 'chr', 'motifstart', 'motifend', 'name', 'score', 'pval', 'strand']):
    
    tissue_cell_allassays = {}
    tissue_cell_assays = {}
    col_list = []
    
    col_list.extend(motif_cols)
    
    with open(tissue_cell_mappings_file, 'r') as ifile:
        lines = ifile.readlines()
        for line in lines:
            if line.startswith('//') or line.startswith('#') or '=' not in line:
                continue
            sl = line.strip().split(key_value_sep)
            key_value = '_'.join(sl[0].strip().replace('(','').replace(')','').replace('-','__').split())
            if key_value not in tissue_cell_assays.keys():
                tissue_cell_assays[key_value] = {}
                tissue_cell_allassays[key_value] = {}
                for assay_name in assay_names:
                    assay_name = '_'.join(assay_name.strip().replace('(','').replace(')','').replace('-','__').split())
                    tissue_cell_allassays[key_value][assay_name] = 'NaN'
                        
            for s in sl[1].split(values_sep):
                cell = s.split(':')[0]
                if cell in cell_assays.keys():
                    cell_updated_name = '_'.join(cell.replace('(','').replace(')','').replace('-','__').split())
                    if ':' in s:
                        assay = s.split(':')[1]
                        if assay in cell_assays[cell]:
                            assay_updated_name = '_'.join(assay.replace('(','').replace(')','').replace('-','__').split())
                            try:
                                tissue_cell_assays[key_value][assay_updated_name].append((cell_updated_name+'___'+assay_updated_name).encode('ascii','ignore').lower())
                            except KeyError:
                                tissue_cell_assays[key_value][assay_updated_name] = [(cell_updated_name+'___'+assay_updated_name).encode('ascii','ignore').lower()]
                            col_list.append((cell_updated_name+'___'+assay_updated_name).encode('ascii','ignore').lower())
                    else:
                        for assay in cell_assays[cell]:
                            assay_updated_name = '_'.join(assay.replace('(','').replace(')','').replace('-','__').split())
                            try:
                                tissue_cell_assays[key_value][assay_updated_name].append((cell_updated_name+'___'+assay_updated_name).encode('ascii','ignore').lower())
                            except KeyError:
                                tissue_cell_assays[key_value][assay_updated_name] = [(cell_updated_name+'___'+assay_updated_name).encode('ascii','ignore').lower()]
                            col_list.append((cell_updated_name+'___'+assay_updated_name).encode('ascii','ignore').lower())
                                
    return col_list, tissue_cell_assays, tissue_cell_allassays

def db_setup_tissues(tissue_cell_assays, 
                     assay_cells_datatypes, 
                     db_name, db_user_name, db_host_name,
                     motif_cols = ['mid INTEGER']):
    
    conn = psycopg2.connect("dbname={} user={} host={}".format(db_name, db_user_name, db_host_name))
    curs = conn.cursor()
    tissue_cols = {}
    for tissue in sorted(tissue_cell_assays.keys()):
        fields = []
        fields.extend(motif_cols)
        
        tissue_cols[tissue] = []
        for c in motif_cols:
            tissue_cols[tissue].append(c.split(' ')[0])
        for assay in sorted(tissue_cell_assays[tissue].keys()):
            fields.append(assay + ' ' + assay_cells_datatypes[assay.replace('__', '-')])
            tissue_cols[tissue].append(assay.encode('ascii','ignore'))
        #create a table for each : tissue
        #print 'select count(mid) from {} where mid<= 45000000;'.format(tissue)
        #curs.execute('select count(mid) from {} where mid<= 45000000;'.format(tissue))
        #print curs.fetchone()[0]
        curs.execute("DROP TABLE IF EXISTS {}".format(tissue))
        create_table_stmt = "CREATE TABLE IF NOT EXISTS {} ({});".format(tissue, ' ,'.join(fields))
        print create_table_stmt
        curs.execute(create_table_stmt)
    curs.close()
    conn.commit()
    conn.close()
    return tissue_cols

def insert_into_tissues(selected_rows, tissue_cell_assays, tissue_cell_allassays, assay_names,
                           cols_to_write_to,
                           cols_to_write_to_allassays,thread_num,
                           db_name, db_user_name, db_host_name):
    
    print "Thread {} has started".format(thread_num)
    conn = open_connection(db_name, db_user_name, db_host_name)
    curs_for_insertion = conn.cursor()
    
    tissues_values = {}
    tissues_fields = {}
    for tissue in sorted(tissue_cell_allassays.keys()):
        tissues_fields[tissue] = ['mid']
        tissues_values[tissue] = []
        for assay in sorted(tissue_cell_allassays[tissue].keys()):
            tissues_fields[tissue].append(assay.encode('ascii','ignore').lower())
        
    #values_to_write = []
    t_process = time.time()
    for row in selected_rows:
        #value_current_row = [row['mid']]
        #reset tissue_cell_allassays for every new row
        for tissue in sorted(tissue_cell_allassays.keys()):
            for assay in sorted(tissue_cell_allassays[tissue].keys()):
                tissue_cell_allassays[tissue][assay]='NaN'
                
        for tissue in sorted(tissue_cell_assays.keys()):
            for assay in sorted(tissue_cell_assays[tissue].keys()):
                values = []
                for col in tissue_cell_assays[tissue][assay]:
                    #col = col.encode('ascii','ignore').lower()
                    #if row[col]!="NaN" and row[col]!=float('nan') and row[col]!='nan':
                    if row[col] not in ["NaN", float('NaN'), float('nan'), 'nan']:
                        try:
                            if not math.isnan(float(row[col])):
                                values.append(float(row[col]))
                        except ValueError:
                            values.append(row[col])              
                
                value = 'NaN'
                if len(values)>0:
                    try:
                        value = sum(values)/float(len(values))
                    except TypeError:
                        value = Counter(values).most_common(1)[0][0]
                    
                tissue_cell_allassays[tissue][assay]=value
                #value_current_row.append(value)
                
        #values_to_write.append(value_current_row)
        #impute missing values
        for assay in assay_names:
            assay = '_'.join(assay.replace('(','').replace(')','').replace('-','__').split())
            tissues_with_NaN_values = []
            tissues_with_values = []
            for tissue in tissue_cell_allassays.keys():
                if tissue_cell_allassays[tissue][assay]=="NaN":
                    tissues_with_NaN_values.append(tissue)
                else:
                    try:
                        tissues_with_values.append(float(tissue_cell_allassays[tissue][assay]))
                    except ValueError:
                        tissues_with_values.append(tissue_cell_allassays[tissue][assay])
            
            if len(tissues_with_NaN_values)>0:
                if len(tissues_with_values)>=4:
                    value = 'NaN'
                    try:
                        value = sum(tissues_with_values)/float(len(tissues_with_values))
                    except TypeError:
                        value = Counter(tissues_with_values).most_common(1)[0][0]
                    if value!='NaN':
                        for tissue in tissues_with_NaN_values:
                            tissue_cell_allassays[tissue][assay] = value
        for tissue in sorted(tissue_cell_allassays.keys()):
            values_selected_row = [row['mid']]
            for assay in sorted(tissue_cell_allassays[tissue].keys()):
                values_selected_row.append(tissue_cell_allassays[tissue][assay])
            tissues_values[tissue].append(values_selected_row)
    print 't_process (func): ', time.time()-t_process
    
    #insert all collected values to their perspective tables/tissues
    t_insert = time.time()
    for tissue in sorted(tissues_values.keys()):
        s_chars = ','.join('%s' for i in range(0, len(tissues_fields[tissue])))
        field_names=','.join(tissues_fields[tissue])
        
        #curs_for_insertion.executemany('insert into {table_name} ({field_names}) values ({values})'.format(table_name='temp'+tissue, field_names=field_names, values=s_chars), tissues_values[tissue])
        '''execute_values(cur = curs_for_insertion, 
                       sql = 'insert into {table_name} ({field_names}) values (%s)'.format(table_name=tissue, field_names=','.join(tissues_fields[tissue])), 
                       argslist = tissues_values[tissue],
                       template = None,
                       page_size = 100)
        '''
        t_tissues_values = tuple(tissues_values[tissue])
        dataText = ','.join('('+curs_for_insertion.mogrify(s_chars, row) + ')' for row in t_tissues_values)
        curs_for_insertion.execute('insert into {table_name} ({field_names}) values {values}'.format(table_name=tissue, field_names=field_names, values=dataText)) 
        
    print 't_insert (func): ', time.time()-t_insert
    conn.commit()
    del tissues_values
    del tissues_fields
    print "Thread {} is done".format(thread_num)
    curs_for_insertion.close()
    conn.close()
    return

def populate_tissue_values(tissue_cell_assays, tissue_cell_allassays, assay_names, col_list, table_from,
                           run_in_parallel_param, 
                           number_processes_to_run_in_parallel,
                           scored_motifs_overlapping_tracks_files,
                           db_name, db_user_name, db_host_name,
                           cols_to_write_to = [], 
                           cols_to_write_to_allassays = [],
                           number_of_rows_to_load = 100000,
                           ):
    
    for tissue in sorted(tissue_cell_assays.keys()):
        for assay in sorted(tissue_cell_assays[tissue].keys()):
            cols_to_write_to.append(tissue+'___'+assay)
    
    for tissue in sorted(tissue_cell_allassays.keys()):
        for assay in sorted(tissue_cell_allassays[tissue].keys()):
            cols_to_write_to_allassays.append(tissue+'___'+assay)
    
    conn = open_connection(db_name, db_user_name, db_host_name)
    curs_for_count = conn.cursor(name = "countcurs", cursor_factory=DictCursor)
    thread_num = 0
    curs_for_count.execute('select count(posrange) from {}'.format(table_from))
    num_rows = int(curs_for_count.fetchone()[0])#85459976
    curs_for_count.close()
    print 'total number of rows to be inserted: ', num_rows
    curs_for_selection = conn.cursor(name = "selectioncurs", cursor_factory=DictCursor)
    curs_for_selection.itersize = number_of_rows_to_load
    t_for_select =  time.time()
    curs_for_selection.execute('select {} from {}'.format(','.join(col_list), table_from))
    print 't_to_select: ', time.time() - t_for_select
    t_for_fetch = time.time()
    selected_rows = curs_for_selection.fetchmany(size=number_of_rows_to_load)
    print 't_to_fetch: ', time.time()-t_for_fetch
    print "Selected {} rows for insertion.".format(str(number_of_rows_to_load))
    t_jobset = time.time()
    if run_in_parallel_param and len(scored_motifs_overlapping_tracks_files)>1:
        print 'Running in parallel'
        i = 0
        num_cores = number_processes_to_run_in_parallel
        p = Pool(number_processes_to_run_in_parallel)
        while i<num_cores:
            p.apply_async(insert_into_tissues, args=(selected_rows, tissue_cell_assays, tissue_cell_allassays, assay_names,
                                   cols_to_write_to, cols_to_write_to_allassays, thread_num,
                                   db_name, db_user_name, db_host_name))
            num_rows -=len(selected_rows)
            print 'num_rows remaining: ', num_rows
            
            selected_rows = []
            selected_rows = curs_for_selection.fetchmany(size=number_of_rows_to_load)
            
            if len(selected_rows)<=0:
                p.close()
                p.join()
                break
            if i==num_cores-1:
                p.close()
                p.join()
                print 't_jobset: ', time.time()-t_jobset
                t_jobset = time.time()
                p = Pool(number_processes_to_run_in_parallel)
                i=0
            i+=1
            thread_num+=1
    else:
        print 'Running sequentially'
        while selected_rows:
            t_to_insert = time.time()
            insert_into_tissues(selected_rows, tissue_cell_assays, tissue_cell_allassays, assay_names,
                                   cols_to_write_to, cols_to_write_to_allassays, thread_num,
                                   db_name, db_user_name, db_host_name)
            
            print 't_to_insert: ', time.time()-t_to_insert
            num_rows -=len(selected_rows)
            print 'num_rows remaining: ', num_rows
            t_for_fetch = time.time()
            selected_rows = []
            selected_rows = curs_for_selection.fetchmany(size=number_of_rows_to_load)
            thread_num+=1
            print 't_to_fetch: ', time.time()-t_for_fetch
            if num_rows<=0:
                print "All rows are processed and inserted from {} into tissue tables".format(table_from)
                sys.exit(0)
    
    if num_rows==0:
        print "All rows are processed and inserted from {} into tissue tables".format(table_from)
    else:
        print "Warning!!! There are {} remaining rows to be inserted".format(num_rows)
    curs_for_selection.close()
    conn.close()
    return
    
def update_table(conn, cell_table, 
                 update_chr_col = True, 
                 chr_col='chr', 
                 update_posrange_col = True, 
                 posrange_col = 'posrange', 
                 motifstart_col='motifstart', 
                 motifend_col='motifend', 
                 pval_score_update=True, 
                 score_col_name='score', 
                 pval_col_name = 'pval'):
    curs = conn.cursor()
    curs.execute("alter table {} add column IF NOT EXISTS {} int4range".format(cell_table, posrange_col))
    
    if update_chr_col:
        curs.execute("update {0} set {1} = replace(replace(replace(replace({1},'chr', ''), 'X', '23'), 'Y', '24'), 'M', '25');".format(cell_table, chr_col))
        curs.execute("alter table {0} alter column {1} SET DATA TYPE integer USING {1}::integer".format(cell_table, chr_col))
    
    if update_posrange_col:
        curs.execute("update {0} set {1} = int4range({2}, {3}, '[]')".format(cell_table, posrange_col, motifstart_col, motifend_col))
        curs.execute("alter table {0} alter column {1} SET NOT NULL".format(cell_table, posrange_col))
    
    if pval_score_update:
        curs.execute("alter table {0} add column IF NOT EXISTS {1} real".format(cell_table, pval_col_name))
        curs.execute("update {0} set {1} = split_part({2}, 'P', 2)::real".format(cell_table, pval_col_name, score_col_name))
        
        curs.execute("update {0} set {1} = replace(split_part({1}, 'P', 1), 'S', '')".format(cell_table,score_col_name))
        curs.execute("alter table {0} alter column {1} SET DATA TYPE real USING {1}::real".format(cell_table, score_col_name))
    conn.commit()
    
    return

def create_index(conn, cell_table, index_name='indexposrange', index_method = 'gist', index_cols = 'posrange'):
    curs = conn.cursor()
    curs.execute("DROP INDEX IF EXISTS {}".format(index_name))
    creat_index_stmt= "CREATE INDEX IF NOT EXISTS {} ON {} using {} ({})".format(index_name, cell_table, index_method, index_cols)
    print creat_index_stmt
    curs.execute(creat_index_stmt)
    conn.commit()
    #conn.close()
    return

def table_contains_data(conn, table_name):
    
    curs = conn.cursor()
    curs.execute('select chr from {} limit 1'.format(table_name))
    if curs.fetchone() is not None:
        print '{} contains data'.format(table_name)
        return True
    else:
        return False

def create_table_stmt_parallel(db_name, db_user_name, db_host_name,
                               tissue, tissuecols, tissuemotifsimputed):
    
    conn = open_connection(db_name, db_user_name, db_host_name)
    curs = conn.cursor()
    print "Loading for", tissue
    print tissuecols
    print "create table if not exists {0} as (select {1} from {2})".format(tissue, tissuecols, tissuemotifsimputed)
    curs.execute("create table if not exists {0} as (select {1} from {2})".format(tissue, tissuecols, tissuemotifsimputed))
    curs.execute('alter table {0} add column if not exists mid serial unique references motifs(mid);'.format(tissue))
    conn.commit()
    print "Created", tissue
    curs.close()
    conn.close()

def populate_table_per_tissue(tissues_table_name, db_name, db_user_name, db_host):
    
    conn = open_connection(db_name, db_user_name, db_host)
    tissues_col_names = get_col_names_from_table(tissues_table_name, conn)
    conn.close()
    tissues = {}
    for col in tissues_col_names:
        if '___' in col:
            if col.split('___')[0] not in tissues.keys():
                tissues[col.split('___')[0]] = [col + ' as ' + col.split('___')[1]]
            else:
                tissues[col.split('___')[0]].append(col + ' as ' + col.split('___')[1])
    p = Pool()
    for tissue in tissues:
        #p.apply_async(create_table_stmt_parallel, args= (db_name, tissue, ','.join(tissues[tissue]), 'tissuemotifsimputed'))
        create_table_stmt_parallel(db_name, db_user_name, db_host,
                                    tissue, ','.join(tissues[tissue]), 'tissuemotifsimputed')   
    p.close()
    p.join()
    print "All tables are added successfully"

def split_motifs_parallel(db_name, db_user_name, db_host, motifs_table, chr, motif_cols):
    conn = open_connection(db_name, db_user_name, db_host)
    curs = conn.cursor()
    new_table_name = "chr"+str(chr)+"motifs"
    print new_table_name
    curs.execute('drop table if exists {0};'.format(new_table_name))
    print 'create table if not exists {0} as (select {1} from {2} where chr={3});'.format(new_table_name, ','.join(motif_cols), motifs_table, chr)
    curs.execute('create table if not exists {0} as (select {1} from {2} where chr={3});'.format(new_table_name, ','.join(motif_cols), motifs_table, chr))
    curs.execute('create index if not exists {0}mid on {1} using btree(mid);'.format("chr"+str(chr), new_table_name))
    curs.execute('create index if not exists {0}posrange on {1} using gist(posrange);'.format("chr"+str(chr), new_table_name))
    conn.commit()
    curs.close()
    conn.close()

def split_motifs_table_by_chr(motifs_table, motif_cols, db_name):
    
    chr_names = range(1,26)
    print db_name, motifs_table
    p = Pool()
    for chr in chr_names:
        p.apply_async(split_motifs_parallel, args = (db_name, motifs_table, chr, motif_cols))
        #split_motifs_parallel(db_name, motifs_table, chr)
    p.close()
    p.join()
    print 'All tables are created'
    return

#to update motif tables
def update_motif_pos_pertable(db_name, db_user_name, db_host_name, chr_table):

    conn = open_connection(db_name, db_user_name, db_host_name)
    curs = conn.cursor()
    cmd =  'update {} set posrange=int4range(lower(posrange)+1,upper(posrange)+1), motifstart=motifstart+1,motifend=motifend+1;'.format(chr_table)
    print cmd
    curs.execute(cmd)
    conn.commit()
    curs.close()
    conn.close()
    
    return

def update_motif_pos(db_name, db_user_name, db_host_name):
    
    p = Pool(8)
    chr_names = range(1,26)
    for chr in chr_names:
        chr_table = "chr"+str(chr)+"motifs"
        p.apply_async(update_motif_pos_pertable, args=(db_name, chr_table))
    p.close()
    p.join()
    print "updated all the tables"
    return


def generate_db(db_name,
                cell_table,
                db_user_name,
                db_host_name,
                cells_assays_dict,
                assay_cells_datatypes,
                cell_assays,
                assay_names,
                tissue_cell_mappings_file,
                run_in_parallel_param,
                number_processes_to_run_in_parallel,
                header,
                scored_motifs_overlapping_tracks_files,
                motif_cols,
                motif_cols_names,
                cell_index_name='indexposrange', cell_index_method = 'gist', cell_index_cols = 'posrange',
                number_of_rows_to_load=50000
        ):
    field_names = db_setup(cells_assays_dict, 
             assay_cells_datatypes, 
             db_name = db_name, 
             cell_table=cell_table, 
             db_user_name=db_user_name, 
             db_host_name=db_host_name,
             motif_cols=motif_cols,
             motif_col_names=motif_cols_names[1:])#the first column mid is auto-incremental
    
    conn = open_connection(db_name, db_user_name, db_host_name)
    
    if not table_contains_data(conn, cell_table): #that is to avoid writing over an already existing table content
        print "Inserting data into: ", cell_table
        insert_into_db(field_names, db_name = db_name, 
                       db_user_name=db_user_name, 
                       db_host_name=db_host_name, 
                       cell_table=cell_table, 
                       scored_motifs_overlapping_tracks_files=scored_motifs_overlapping_tracks_files, 
                       header=header,
                       run_in_parallel_param=run_in_parallel_param,
                       number_processes_to_run_in_parallel=number_processes_to_run_in_parallel)# dir_to_import=params['motifs_overlapping_tracks_output_dir'], keyword_to_check="_scored.bed10", header=header)
        print "Creating index on: ", cell_table
        create_index(conn, cell_table, index_name=cell_index_name, index_method = cell_index_method, index_cols = cell_index_cols)
    close_connection(conn)
    
    process_tissues = True
    #write results to the tissues (based on cell motifs) table
    if process_tissues:
        print 'Creating tissues tables'
        col_list, tissue_cell_assays, tissue_cell_allassays = get_tissue_cell_mappings(cell_assays, 
                                                                                       assay_names, 
                                                                                       tissue_cell_mappings_file,  
                                                                                       motif_cols = motif_cols_names)
        tissue_cols = db_setup_tissues(tissue_cell_allassays, 
                                       assay_cells_datatypes, 
                                       motif_cols = ['mid INTEGER'], 
                                       db_name = db_name, db_user_name=db_user_name, db_host_name=db_host_name)
        col_list.append('mid')
        print 'Inserting data into tissues tables'
        populate_tissue_values(tissue_cell_assays, 
                               tissue_cell_allassays, 
                               assay_names, col_list, table_from=cell_table, 
                               run_in_parallel_param=run_in_parallel_param, 
                               number_processes_to_run_in_parallel=number_processes_to_run_in_parallel, 
                               scored_motifs_overlapping_tracks_files=scored_motifs_overlapping_tracks_files, 
                               db_name = db_name, db_user_name=db_user_name, db_host_name=db_host_name,
                               number_of_rows_to_load = number_of_rows_to_load,
                               )
        
        print "Creating index on tissues tables"
        conn = open_connection(db_name, db_user_name, db_host_name)
        for tissue_table in sorted(tissue_cols.keys()):
            
            create_index(conn, tissue_table, index_name='index'+tissue_table+'mid', index_method = 'btree', index_cols = 'mid')
        conn.close()
    
    #split motif table per chr
    split_motifs = False
    if split_motifs:
        split_motifs_table_by_chr(motifs_table=cell_table, motif_cols=motif_cols_names, db_name=db_name)
    
    update_motif_positions = False
    if update_motif_positions:
        update_motif_pos(db_name, db_user_name, db_host_name)
        