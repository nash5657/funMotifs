'''
Created on 6 Oct 2017

@author: husensofteng
'''

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
                     motif_cols = ['mid INTEGER'],
                     tissues_fscores_table="tissues_fscores_table",
                     motif_cols_for_tissues_fscores_table = ['mid INTEGER']):
    
    conn = psycopg2.connect("dbname={} user={} host={}".format(db_name, db_user_name, db_host_name))
    curs = conn.cursor()
    tissue_cols = {}
    fields_tissues = []
    fields_tissues.extend(motif_cols_for_tissues_fscores_table)
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
        
        fields_tissues.append(tissue + " numeric")#for the fscores table of all tissues
    
    curs.execute("DROP TABLE IF EXISTS {}".format(tissues_fscores_table))
    create_table_stmt = "CREATE TABLE IF NOT EXISTS {} ({});".format(tissues_fscores_table, ' ,'.join(fields_tissues))
    print create_table_stmt
    curs.execute(create_table_stmt)
    curs.close()
    conn.commit()
    conn.close()
    return tissue_cols

def get_score_from_value(value, assay, feature_weights_dict):
    score = 0.0
    try:
        score  = feature_weights_dict[value]#where the value/label is present in the weights conf file
    except KeyError:
        try:
            if float(value)>0:
                score = feature_weights_dict[assay]#where the assay name is present in the weights conf file
        except ValueError:#for number of items where each item adds a unit of the weight
            if value!="" and value!=" ":
                try:
                    score = feature_weights_dict[assay]*len(value.split(','))
                except KeyError:
                    return score
        except KeyError:
            return score 
    return score

def insert_into_tissues(selected_rows, tissue_cell_assays, tissue_cell_allassays, assay_names,
                           cols_to_write_to,
                           cols_to_write_to_allassays,thread_num, 
                           feature_weights_dict, 
                           db_name, db_user_name, db_host_name, tissues_fscores_table):
            
    print "Thread {} has started".format(thread_num)
    conn = DBUtilities.open_connection(db_name, db_user_name, db_host_name)
    curs_for_insertion = conn.cursor()
    
    tissues_values = {}
    tissues_fields = {}
    for tissue in sorted(tissue_cell_allassays.keys()):
        tissues_fields[tissue] = ['mid', 'fscore']
        tissues_values[tissue] = []
        for assay in sorted(tissue_cell_allassays[tissue].keys()):
            tissues_fields[tissue].append(assay.encode('ascii','ignore').lower())
    fscores_per_tissues_allrows = []
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
        
        fscores_per_tissues = [row['mid']]
        for tissue in sorted(tissue_cell_allassays.keys()):
            values_selected_row = [row['mid'], 0.0]
            fscore = 0.0
            for assay in sorted(tissue_cell_allassays[tissue].keys()):
                values_selected_row.append(tissue_cell_allassays[tissue][assay])
                
                #compute the final score
                value = tissue_cell_allassays[tissue][assay]
                fscore+=get_score_from_value(value, assay, feature_weights_dict)
                
            values_selected_row[1]=fscore
            tissues_values[tissue].append(values_selected_row)
            #for the tissues_fscores table
            fscores_per_tissues.append(fscore)
        fscores_per_tissues_allrows.append(fscores_per_tissues)
    print 't_process (func): ', time.time()-t_process
    
    #insert all collected values to their respective tables/tissues
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
        
    #insert into tissues_fscores table
    tissues_names_for_fscores = ['mid']
    tissues_names_for_fscores.extend(sorted(tissue_cell_allassays.keys()))
    print ', '.join(tissues_names_for_fscores)
    s_chars_for_fscores = ','.join('%s' for i in range(0, len(tissues_names_for_fscores)))
    print s_chars_for_fscores
    print fscores_per_tissues_allrows
    t_tissues_fscores_values = tuple(fscores_per_tissues_allrows)
    fscores_per_tissues_dataText = ','.join('('+curs_for_insertion.mogrify(s_chars_for_fscores, row) + ')' for row in t_tissues_fscores_values)
    print fscores_per_tissues_dataText
    curs_for_insertion.execute('insert into {table_name} ({field_names}) values {values}'.format(table_name=tissues_fscores_table, field_names=', '.join(tissues_names_for_fscores), values=fscores_per_tissues_dataText)) 
        
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
                           feature_weights_dict,
                           db_name, db_user_name, db_host_name,
                           tissues_fscores_table,
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
                                   cols_to_write_to, cols_to_write_to_allassays, thread_num, feature_weights_dict,
                                   db_name, db_user_name, db_host_name, tissues_fscores_table))
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
                                   cols_to_write_to, cols_to_write_to_allassays, thread_num, feature_weights_dict,
                                   db_name, db_user_name, db_host_name, tissues_fscores_table)
            
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
                break
    
    if num_rows==0:
        print "All rows are processed and inserted from {} into tissue tables".format(table_from)
    else:
        print "Warning!!! There are {} remaining rows to be inserted".format(num_rows)
    curs_for_selection.close()
    conn.close()
    return
    
def get_weights_per_feature(annotation_weights_inputfile, skip_negative_weights):
    
    feature_weights = {}
    with open(annotation_weights_inputfile, 'r') as ifile:
        lines = ifile.readlines()
    for l in lines:
        if not l.startswith('#') and not l.startswith('//'):
            sl = l.split('#')[0].split('=')
            if len(sl)==2:
                try:
                    for label in sl[0].strip().split(','):
                        if skip_negative_weights:
                            if float(sl[1].strip())>=0:
                                feature_weights[label.upper()]=float(sl[1].strip())
                        else:
                            feature_weights[label.upper()]=float(sl[1].strip())
                except ValueError:
                    print "no weights is used for sl[0] because the assigned values is not a number"
                    continue
    return feature_weights

def generate_tissue_tables(db_name,
                    cell_table,
                    db_user_name,
                    db_host_name,
                    tissues_fscores_table,
                    assay_cells_datatypes,
                    cell_assays,
                    assay_names,
                    tissue_cell_mappings_file,
                    run_in_parallel_param,
                    number_processes_to_run_in_parallel,
                    scored_motifs_overlapping_tracks_files,
                    motif_cols_names,
                    number_of_rows_to_load,
                    annotation_weights_inputfile,
                    skip_negative_weights,
                    ):
    col_list, tissue_cell_assays, tissue_cell_allassays = get_tissue_cell_mappings(cell_assays, 
                                                                                       assay_names, 
                                                                                       tissue_cell_mappings_file,  
                                                                                       motif_cols = motif_cols_names)
    tissue_cols = db_setup_tissues(tissue_cell_allassays, 
                                   assay_cells_datatypes, 
                                   motif_cols = ['mid INTEGER', 'fscore NUMERIC'], 
                                   db_name = db_name, db_user_name=db_user_name, db_host_name=db_host_name,
                                   tissues_fscores_table=tissues_fscores_table,
                                   motif_cols_for_tissues_fscores_table = ['mid INTEGER'])
    col_list.append('mid')
    feature_weights_dict = get_weights_per_feature(annotation_weights_inputfile=annotation_weights_inputfile,skip_negative_weights=skip_negative_weights)
    print 'Inserting data into tissues tables'
    populate_tissue_values(tissue_cell_assays, 
                           tissue_cell_allassays, 
                           assay_names, col_list, table_from=cell_table, 
                           run_in_parallel_param=run_in_parallel_param, 
                           number_processes_to_run_in_parallel=number_processes_to_run_in_parallel, 
                           scored_motifs_overlapping_tracks_files=scored_motifs_overlapping_tracks_files, 
                           feature_weights_dict=feature_weights_dict,
                           db_name = db_name, db_user_name=db_user_name, db_host_name=db_host_name,
                           tissues_fscores_table=tissues_fscores_table,
                           number_of_rows_to_load = number_of_rows_to_load,
                           )
    
    print "Creating index on tissues tables"
    for tissue_table in sorted(tissue_cols.keys()):
        DBUtilities.create_index(db_name, db_user_name, db_host_name, 
                                 tissue_table, index_name='index'+tissue_table+'mid', index_method = 'btree', index_cols = 'mid')
    DBUtilities.create_index(db_name, db_user_name, db_host_name, 
                                 tissue_table=tissues_fscores_table, index_name='index'+tissues_fscores_table+'mid', index_method = 'btree', index_cols = 'mid')