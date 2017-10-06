'''
Created on 6 Oct 2017

@author: husensofteng
'''
import sys
from collections import Counter
from multiprocessing import Pool
import time
import math
import psycopg2
from psycopg2.extras import DictCursor
import DBUtilities

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
    conn = DBUtilities.open_connection(db_name, db_user_name, db_host_name)
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
    
    conn = DBUtilities.open_connection(db_name, db_user_name, db_host_name)
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
    
def generate_tissue_tables(db_name,
                    cell_table,
                    db_user_name,
                    db_host_name,
                    assay_cells_datatypes,
                    cell_assays,
                    assay_names,
                    tissue_cell_mappings_file,
                    run_in_parallel_param,
                    number_processes_to_run_in_parallel,
                    scored_motifs_overlapping_tracks_files,
                    motif_cols_names,
                    number_of_rows_to_load
                    ):
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
    conn = DBUtilities.open_connection(db_name, db_user_name, db_host_name)
    for tissue_table in sorted(tissue_cols.keys()):
        DBUtilities.create_index(conn, tissue_table, index_name='index'+tissue_table+'mid', index_method = 'btree', index_cols = 'mid')
    conn.close()
