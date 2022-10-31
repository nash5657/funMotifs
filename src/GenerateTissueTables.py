'''
Created on 6 Oct 2017

@author: husensofteng
'''

from collections import Counter
from multiprocessing import Pool
import pathos.multiprocessing as mp
import time
import math
import psycopg2
from psycopg2.extras import DictCursor
import DBUtilities
import pandas as pd
import os
import glob
import json

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
            #print(sl)
            key_value = '_'.join(sl[0].strip().replace('(','').replace(')','').replace('-','__').replace('.','').split())
            print(key_value)
            if key_value not in list(tissue_cell_assays.keys()):
                    tissue_cell_assays[key_value] = {}
                    tissue_cell_allassays[key_value] = {}
                    for assay_name in assay_names:
                        assay_name = '_'.join(assay_name.strip().replace('(','').replace(')','').replace('-','__').split())
                        tissue_cell_allassays[key_value][assay_name] = 'NaN'
            for s in sl[1].split(values_sep):
                    cell = s.split(':')[0]
                    #print(cell)
                    if cell in list(cell_assays.keys()):
                        cell_updated=cell
                        
                        if cell[0].isdigit():
                            #if cell=='22Rv1' or cell=='8988T':
                                cell_updated='a'+cell
                                
                            #if cell=="Ammon's horn":
                            #    cell="Ammons horn"
                            #if cell=="Peyer's patch":
                            #    cell="Peyers patch"

                        #if cell=='22Rv1' or cell=='8988T':
                        #    cell_updated='a'+cell       
                        #if cell=="Ammon's horn":
                        #    cell_updated="Ammons horn"
                        #if cell=="Peyer's patch":
                        #    cell_updated="Peyers patch"
                        cell_updated_name = '_'.join(cell_updated.replace('(','').replace(')','').replace('-','__').replace('.','').replace("'","").split())
                        print(cell_updated_name)
                        print((cell_assays[cell]))
                        if ':' in s:
                            assay = s.split(':')[1]
         
                            if assay in cell_assays[cell]:
                                assay_updated_name = '_'.join(assay.replace('(','').replace(')','').replace('-','__').replace('.','').split())
                                try:
                                    tissue_cell_assays[key_value][assay_updated_name].append((cell_updated_name+'___'+assay_updated_name).lower())
                                except KeyError:
                                    tissue_cell_assays[key_value][assay_updated_name] = [(cell_updated_name+'___'+assay_updated_name).lower()]
                                col_list.append((cell_updated_name+'___'+assay_updated_name))
                        else:
                            for assay in cell_assays[cell]:
                                assay_updated_name = '_'.join(assay.replace('(','').replace(')','').replace('-','__').replace('.','').split())
                                try:
                                    tissue_cell_assays[key_value][assay_updated_name].append((cell_updated_name+'___'+assay_updated_name).lower())
                                except KeyError:
                                    tissue_cell_assays[key_value][assay_updated_name] = [(cell_updated_name+'___'+assay_updated_name).lower()]
                                col_list.append((cell_updated_name+'___'+assay_updated_name))

                                
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
            tissue_cols[tissue].append(assay)
        #create a table for each : tissue
        #print 'select count(mid) from {} where mid<= 45000000;'.format(tissue)
        #curs.execute('select count(mid) from {} where mid<= 45000000;'.format(tissue))
        #print curs.fetchone()[0]
        curs.execute("DROP TABLE IF EXISTS {}".format(tissue))
        create_table_stmt = "CREATE TABLE IF NOT EXISTS {} ({});".format(tissue, ' ,'.join(fields))
        curs.execute(create_table_stmt)
        
        fields_tissues.append(tissue + " numeric")#for the fscores table of all tissues
    
    curs.execute("DROP TABLE IF EXISTS {}".format(tissues_fscores_table))
    create_table_stmt = "CREATE TABLE IF NOT EXISTS {} ({});".format(tissues_fscores_table, ' ,'.join(fields_tissues))
    curs.execute(create_table_stmt)
    curs.close()
    conn.commit()
    conn.close()
    return tissue_cols

def get_score_from_value(value, assay, feature_weights_dict):
    score = 0.0
    try:
        score  = feature_weights_dict[value.upper()]#where the value/label is present in the weights conf file
    except (KeyError,AttributeError):
        try:
            if float(value)>0.0:
                if assay.lower()!="numothertfbinding":
                    score = feature_weights_dict[assay.upper()]#where the assay name is present in the weights conf file
                else:
                    v = float(value)
                    if v>3.0:
                        v = 3.0
                    score = feature_weights_dict[assay.upper()]*v
        except (KeyError, ValueError):
            return score 
    return score

def insert_into_tissues(selected_rows, tissue_cell_assays, tissue_cell_allassays, assay_names,
                           cols_to_write_to,
                           cols_to_write_to_allassays,thread_num, 
                           feature_weights_dict, 
                           db_name, db_user_name, db_host_name, tissues_fscores_table):
    
    print(("Thread {} has started".format(thread_num)))
    conn = DBUtilities.open_connection(db_name, db_user_name, db_host_name)
    curs_for_insertion = conn.cursor()
    
    tissues_values = {}
    tissues_fields = {}
    for tissue in sorted(tissue_cell_allassays.keys()):
        tissues_fields[tissue] = ['mid', 'fscore']
        tissues_values[tissue] = []
        for assay in sorted(tissue_cell_allassays[tissue].keys()):
            tissues_fields[tissue].append(assay.lower())
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
                    # TODO: check if statement below and write else statement
                    if col in row.keys():
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
            for tissue in list(tissue_cell_allassays.keys()):
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
                    try:#FIX: For peaks (ChIP-seq & DNase-seq) don't take the average, rather check if there are more than 4 samples that have a peak then set it to 1 otherwise set it to zero, currently if only one cell has a peak with value x and all other cells have no peak for a motif the imputed value would still be larger than zero. 
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
    print(('t_process (func): ', time.time()-t_process))
    
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
        dataText = ','.join('('+curs_for_insertion.mogrify(s_chars, row).decode("utf-8") + ')' for row in t_tissues_values)
        curs_for_insertion.execute('insert into {table_name} ({field_names}) values {values}'.format(table_name=tissue, field_names=field_names, values=dataText)) 
        
    #insert into tissues_fscores table
    tissues_names_for_fscores = ['mid']
    tissues_names_for_fscores.extend(sorted(tissue_cell_allassays.keys()))
    s_chars_for_fscores = ','.join('%s' for i in range(0, len(tissues_names_for_fscores)))
    t_tissues_fscores_values = tuple(fscores_per_tissues_allrows)
    fscores_per_tissues_dataText = ','.join('('+curs_for_insertion.mogrify(s_chars_for_fscores, row).decode("utf-8") + ')' for row in t_tissues_fscores_values)
    curs_for_insertion.execute('insert into {table_name} ({field_names}) values {values}'.format(table_name=tissues_fscores_table, field_names=', '.join(tissues_names_for_fscores), values=fscores_per_tissues_dataText)) 
    
    print(('t_insert (func): ', time.time()-t_insert))
    conn.commit()
    del tissues_values
    del tissues_fields
    print(("Thread {} is done".format(thread_num)))
    curs_for_insertion.close()
    conn.close()
    return


def insert_into_tissues_from_file(scored_motifs_overlapping_tracks_files_tissue, header_scored_lower, mid_value, tissue_cell_assays_file, tissue_cell_allassays_file, assay_names,
                           cols_to_write_to,
                           cols_to_write_to_allassays,thread_num, 
                           feature_weights_dict_file, 
                           db_name, db_user_name, db_host_name, tissues_fscores_table):
    
    
    print(("Thread {} has started".format(thread_num)))
    
    
    with open(tissue_cell_assays_file, 'r') as tissue_cell_assays_infile:
        tissue_cell_assays = json.load(tissue_cell_assays_infile)
        
    with open(tissue_cell_allassays_file, 'r') as tissue_cell_allassays_infile:
        tissue_cell_allassays = json.load(tissue_cell_allassays_infile)
    
    with open(feature_weights_dict_file, 'r') as feature_weights_dict_infile:
        feature_weights_dict = json.load(feature_weights_dict_infile)
        
    #with open(assay_names_file, 'r') as assay_names_infile:
    #    assay_names = assay_names_infile.readline().strip()
    #    return json.loads(assay_names)
        

        
        
    conn = DBUtilities.open_connection(db_name, db_user_name, db_host_name)
    curs_for_insertion = conn.cursor()
    
    tissues_values = {}
    tissues_fields = {}
    for tissue in sorted(tissue_cell_allassays.keys()):
        print(tissue)
        tissues_fields[tissue] = ['mid', 'fscore']
        tissues_values[tissue] = []
        for assay in sorted(tissue_cell_allassays[tissue].keys()):
            tissues_fields[tissue].append(assay.lower())
    
    fscores_per_tissues_allrows = []
    #values_to_write = []
    t_process = time.time()
    with open(scored_motifs_overlapping_tracks_files_tissue) as data_infile:
        l = data_infile.readline()
        i=0
        while l:
            row = l.strip().split('\t')
            mid = mid_value[i]
            print(('mid' + str(mid)))
            #print(row)
    #for row in selected_rows:
        #value_current_row = [row['mid']]
        #reset tissue_cell_allassays for every new row
            for tissue in sorted(tissue_cell_allassays.keys()):
                for assay in sorted(tissue_cell_allassays[tissue].keys()):
                    tissue_cell_allassays[tissue][assay]='NaN'
                    
            for tissue in sorted(tissue_cell_assays.keys()):
                for assay in sorted(tissue_cell_assays[tissue].keys()):
                    values = []
                    for assay_cell in tissue_cell_assays[tissue][assay]:
                        #col = col.encode('ascii','ignore').lower()
                        #if row[col]!="NaN" and row[col]!=float('nan') and row[col]!='nan':
                        assay_col_id = header_scored_lower.index(assay_cell)
                        #print(col_id)
                        row_assay_col = row[assay_col_id]
                        
                        if row_assay_col not in ["NaN", float('NaN'), float('nan'), 'nan']:
                            #print(row_assay_col)
                            try:
                                if not math.isnan(float(row_assay_col)):
                                    values.append(float(row_assay_col))
                            except ValueError:
                                values.append(row_assay_col)              
                    
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
                for tissue in list(tissue_cell_allassays.keys()):
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
                        try:#FIX: For peaks (ChIP-seq & DNase-seq) don't take the average, rather check if there are more than 4 samples that have a peak then set it to 1 otherwise set it to zero, currently if only one cell has a peak with value x and all other cells have no peak for a motif the imputed value would still be larger than zero. 
                            value = sum(tissues_with_values)/float(len(tissues_with_values))
                        except TypeError:
                            value = Counter(tissues_with_values).most_common(1)[0][0]
                        if value!='NaN':
                            for tissue in tissues_with_NaN_values:
                                tissue_cell_allassays[tissue][assay] = value
            
            #fscores_per_tissues = [row['mid']]
            fscores_per_tissues = [mid]
            for tissue in sorted(tissue_cell_allassays.keys()):
                #values_selected_row = [row['mid'], 0.0]
                values_selected_row = [mid, 0.0]
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
            print(fscores_per_tissues)
            fscores_per_tissues_allrows.append(fscores_per_tissues)
            i+=1 
            l = data_infile.readline()
    print(('t_process (func): ', time.time()-t_process))
    
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
        # BUG: here is the bug. "can only concatenate str (not "bytes") to str"
        dataText = ','.join('('+curs_for_insertion.mogrify(s_chars, row) + ')' for row in t_tissues_values)
        try:
            curs_for_insertion.execute('insert into {table_name} ({field_names}) values ({values})'.format(table_name=tissue, field_names=field_names, values=dataText)) 
            tissues_names_for_fscores = ['mid']
            tissues_names_for_fscores.extend(sorted(tissue_cell_allassays.keys()))
            s_chars_for_fscores = ','.join('%s' for i in range(0, len(tissues_names_for_fscores)))
            t_tissues_fscores_values = tuple(fscores_per_tissues_allrows)
            fscores_per_tissues_dataText = ','.join('('+curs_for_insertion.mogrify(s_chars_for_fscores, row) + ')' for row in t_tissues_fscores_values)
            curs_for_insertion.execute('insert into {table_name} ({field_names}) values {values}'.format(table_name=tissues_fscores_table, field_names=', '.join(tissues_names_for_fscores), values=fscores_per_tissues_dataText)) 
    
            print(('t_insert (func): ', time.time()-t_insert))
            conn.commit()
            print("Data inserted into DB")
        except:
            print("Thread coundn't insert to DB ")
    #insert into tissues_fscores table
    
    del tissues_values
    del tissues_fields
    print(("Thread {} is done".format(thread_num)))
    
    curs_for_insertion.close()
    conn.close()
    return

def insert_into_tissues_per_tissue(scored_motifs_overlapping_tracks_files, tissue_cell_assays, tissue_cell_allassays, assay_names,
                           cols_to_write_to,
                           cols_to_write_to_allassays,thread_num, 
                           feature_weights_dict, 
                           db_name, db_user_name, db_host_name, tissues_fscores_table):
    
    print(("Thread {} has started".format(thread_num)))
    conn = DBUtilities.open_connection(db_name, db_user_name, db_host_name)
    curs_for_insertion = conn.cursor()
    
    tissues_values = {}
    tissues_fields = {}
    tissues_fields[tissue] = ['mid', 'fscore']
    tissues_values[tissue] = []
    for assay in sorted(tissue_cell_allassays[tissue].keys()):
        tissues_fields[tissue].append(assay.lower())
    fscores_per_tissues_allrows = []
    #values_to_write = []
    t_process = time.time()
    with open(scored_motifs_overlapping_tracks_files_tissue) as data_infile:
        l = data_infile.readline()
        while l:
            row = l.strip().split('\t')
            for assay in sorted(tissue_cell_allassays[tissue].keys()):
                tissue_cell_allassays[tissue][assay]='NaN'
                
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
                
            l = data_infile.readline()
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
            for tissue in list(tissue_cell_allassays.keys()):
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
                    try:#FIX: For peaks (ChIP-seq & DNase-seq) don't take the average, rather check if there are more than 4 samples that have a peak then set it to 1 otherwise set it to zero, currently if only one cell has a peak with value x and all other cells have no peak for a motif the imputed value would still be larger than zero. 
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
    print(('t_process (func): ', time.time()-t_process))
    
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
    s_chars_for_fscores = ','.join('%s' for i in range(0, len(tissues_names_for_fscores)))
    t_tissues_fscores_values = tuple(fscores_per_tissues_allrows)
    fscores_per_tissues_dataText = ','.join('('+curs_for_insertion.mogrify(s_chars_for_fscores, row) + ')' for row in t_tissues_fscores_values)
    curs_for_insertion.execute('insert into {table_name} ({field_names}) values {values}'.format(table_name=tissues_fscores_table, field_names=', '.join(tissues_names_for_fscores), values=fscores_per_tissues_dataText)) 
    
    print(('t_insert (func): ', time.time()-t_insert))
    conn.commit()
    del tissues_values
    del tissues_fields
    print(("Thread {} is done".format(thread_num)))
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
                           number_of_rows_to_load = 10,
                           ):
    
    for tissue in sorted(tissue_cell_assays.keys()):
        for assay in sorted(tissue_cell_assays[tissue].keys()):
            cols_to_write_to.append(tissue+'___'+assay)
    
    for tissue in sorted(tissue_cell_allassays.keys()):
        for assay in sorted(tissue_cell_allassays[tissue].keys()):
            cols_to_write_to_allassays.append(tissue+'___'+assay)
    
    col_list_lower = [x.lower() for x in col_list] 
    conn = DBUtilities.open_connection(db_name, db_user_name, db_host_name)
    curs_for_count = conn.cursor(name = "countcurs", cursor_factory=DictCursor)
    thread_num = 0
    curs_for_count.execute('select count(posrange) from {}'.format(table_from))
    num_rows = int(curs_for_count.fetchone()[0])#85459976
    curs_for_count.close()
    print(('total number of rows to be inserted: ', num_rows))
    curs_for_rownames = conn.cursor(name = "rownamecurs", cursor_factory=DictCursor)
    curs_for_rownames.execute('SELECT * from {}'.format(table_from))
    first_row = curs_for_rownames.fetchone()
    colnames = [desc[0] for desc in curs_for_rownames.description]
    curs_for_selection = conn.cursor(name = "selectioncurs", cursor_factory=DictCursor)
    curs_for_selection.itersize = number_of_rows_to_load
    t_for_select =  time.time()
    # TODO: check and integrate following line
    col_list_lower = colnames
    curs_for_selection.execute('select {} from {}'.format(','.join(col_list_lower), table_from))
    print(('t_to_select: ', time.time() - t_for_select))
    t_for_fetch = time.time()
    selected_rows = curs_for_selection.fetchmany(size=number_of_rows_to_load)
    print(('t_to_fetch: ', time.time()-t_for_fetch))
    print(("Selected {} rows for insertion.".format(str(number_of_rows_to_load))))
    t_jobset = time.time()
    if run_in_parallel_param and len(scored_motifs_overlapping_tracks_files)>1:
        print('Running in parallel')
        i = 0
        num_cores = number_processes_to_run_in_parallel
        p = Pool(number_processes_to_run_in_parallel)
        while i<num_cores:
            p.apply_async(insert_into_tissues, args=(selected_rows, tissue_cell_assays, tissue_cell_allassays, assay_names,
                                   cols_to_write_to, cols_to_write_to_allassays, thread_num, feature_weights_dict,
                                   db_name, db_user_name, db_host_name, tissues_fscores_table))
            num_rows -=len(selected_rows)
            print(('num_rows remaining: ', num_rows))
            
            selected_rows = []
            selected_rows = curs_for_selection.fetchmany(size=number_of_rows_to_load)
            
            if len(selected_rows)<=0:
                p.close()
                p.join()
                break
            if i==num_cores-1:
                p.close()
                p.join()
                print(('t_jobset: ', time.time()-t_jobset))
                t_jobset = time.time()
                p = Pool(number_processes_to_run_in_parallel)
                i=0
            i+=1
            thread_num+=1
    else:
        print('Running sequentially')
        while selected_rows:
            t_to_insert = time.time()
            insert_into_tissues(selected_rows, tissue_cell_assays, tissue_cell_allassays, assay_names,
                                   cols_to_write_to, cols_to_write_to_allassays, thread_num, feature_weights_dict,
                                   db_name, db_user_name, db_host_name, tissues_fscores_table)
            
            print(('t_to_insert: ', time.time()-t_to_insert))
            num_rows -=len(selected_rows)
            print(('num_rows remaining: ', num_rows))
            t_for_fetch = time.time()
            selected_rows = []
            selected_rows = curs_for_selection.fetchmany(size=number_of_rows_to_load)
            thread_num+=1
            print(('t_to_fetch: ', time.time()-t_for_fetch))
            if num_rows<=0:
                print(("All rows are processed and inserted from {} into tissue tables".format(table_from)))
                break
    
    if num_rows==0:
        print(("All rows are processed and inserted from {} into tissue tables".format(table_from)))
    else:
        print(("Warning!!! There are {} remaining rows to be inserted".format(num_rows)))
    curs_for_selection.close()
    conn.close()
    return
    
    
    
def populate_tissue_values_from_scored_files(tissue_cell_assays, tissue_cell_allassays, assay_names, col_list, 
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
    
    dir_in = os.path.dirname(scored_motifs_overlapping_tracks_files[0])
    tissue_cell_assays_file = dir_in+"/tissue_cell_assays_dict"
    tissue_cell_allassays_file = dir_in+"/tissue_cell_allassays_dict"
    feature_weights_dict_file = dir_in+ "/feature_weights_dict"
    #assay_names_file = dir_in+ "/assay_names"
    
    with open(tissue_cell_assays_file, 'w') as tissue_cell_assays_outfile:
        json.dump(tissue_cell_assays, tissue_cell_assays_outfile)
    
    
    for tissue in sorted(tissue_cell_allassays.keys()):
        for assay in sorted(tissue_cell_allassays[tissue].keys()):
            cols_to_write_to_allassays.append(tissue+'___'+assay)
            
            
    with open(tissue_cell_allassays_file, 'w') as tissue_cell_allassays_outfile:
        json.dump(tissue_cell_allassays, tissue_cell_allassays_outfile)        
        
    with open(feature_weights_dict_file, 'w') as feature_weights_dict_outfile:
        json.dump(feature_weights_dict, feature_weights_dict_outfile) 
        
   # with open(assay_names_file, 'w') as assay_names_outfile:
    #    json.dump(assay_names, assay_names_outfile) 
    
    thread_num = 0
    
    col_list_lower = [x.lower() for x in col_list] 
#     conn = DBUtilities.open_connection(db_name, db_user_name, db_host_name)
#     curs_for_count = conn.cursor(name = "countcurs", cursor_factory=DictCursor)
#     thread_num = 0
#     curs_for_count.execute('select count(posrange) from {}'.format(table_from))
#     num_rows = int(curs_for_count.fetchone()[0])#85459976
#     curs_for_count.close()
#     print('total number of rows to be inserted: ', num_rows)
#     curs_for_selection = conn.cursor(name = "selectioncurs", cursor_factory=DictCursor)
#     curs_for_selection.itersize = number_of_rows_to_load
#     t_for_select =  time.time()
#     curs_for_selection.execute('select {} from {}'.format(','.join(col_list), table_from))
    #print('t_to_select: ', time.time() - t_for_select)
    #t_for_fetch = time.time()
    #selected_rows = curs_for_selection.fetchmany(size=number_of_rows_to_load)
    #print('t_to_fetch: ', time.time()-t_for_fetch)
    #print("Selected {} rows for insertion.".format(str(number_of_rows_to_load)))
    t_jobset = time.time()
    if run_in_parallel_param and len(scored_motifs_overlapping_tracks_files)>1:
        print('Running in parallel')
        i = 0
        #n_lines = 1
       # n_lines_end=number_of_rows_to_load
        num_cores = number_processes_to_run_in_parallel
        
        for file_in in scored_motifs_overlapping_tracks_files:
            #i = 0
            #file_id=0
            #read header of scored motifs file
            with open(file_in) as f:
                header_scored = f.readline().strip('"').split('\t')
            
            header_scored_lower = [x.lower().replace('"','') for x in header_scored] 
            
            #divide files into subfiles
            comm_divide_files ="cat {} | tail -n +2 | split -l {} - {}".format(file_in, str(number_of_rows_to_load), file_in+"_part")
            #print(comm_divide_files)
            os.system(comm_divide_files)
            part_file_list = glob.glob(file_in+'_part*', recursive=False)
            #n_part_file_list = len(part_file_list)
            #print(n_part_file_list)
            
            
            #n_lines = 1
            #n_lines_end=number_of_rows_to_load
            num_cores = number_processes_to_run_in_parallel
            p = mp.Pool(number_processes_to_run_in_parallel)
            #selected_rows_df= pd.read_csv(file_in, nrows=number_of_rows_to_load, sep='\t', dtype=str)
            #selected_rows_df_order = selected_rows_df.reindex(columns=col_list[1::])
            #to list of tuples
            #records = selected_rows_df_order.to_records(index=False)
            #selected_rows = list(records)
            
            #while i<num_cores:
            for part_file in part_file_list: 
                
                part_file_num_lines = sum(1 for _ in open(part_file))
                mid_value = list(range(i,i+part_file_num_lines))
                print(mid_value)
                i = i + part_file_num_lines
                
                res = p.apply_async(insert_into_tissues_from_file, args=[part_file, header_scored_lower, mid_value, tissue_cell_assays_file, tissue_cell_allassays_file, list(assay_names),
                                       cols_to_write_to, cols_to_write_to_allassays, thread_num, feature_weights_dict_file,
                                       db_name, db_user_name, db_host_name, tissues_fscores_table])
                
                
                print((res.get()))
                thread_num+=1
                #num_rows -=len(selected_rows)
                #print('num_rows: ', n_lines_end)
                
               # selected_rows = []
               # selected_rows_df= pd.read_csv(file_in,skiprows=range(n_lines,n_lines_end), nrows=number_of_rows_to_load, sep='\t', dtype=str)
               # selected_rows_df_order = selected_rows_df.reindex(columns=col_list[1::])
                #to list of tuples
                #records = selected_rows_df_order.to_records(index=False)
                #selected_rows = list(records)
                #n_lines_end = n_lines_end + number_of_rows_to_load
                #selected_rows = curs_for_selection.fetchmany(size=number_of_rows_to_load)
                
                #file_id+=1
                #if n_file_list==file_id:
            p.close()
            p.join()
                #   break
                #if i==num_cores-1:
                #    p.close()
                #    p.join()
                #    print('t_jobset: ', time.time()-t_jobset)
                #    t_jobset = time.time()
                #    p = Pool(number_processes_to_run_in_parallel)
                #    i=0
                #i+=1
                
            
            #for file_in_part in  file_list:
             #   os.remove(file_in_part)   
                
    else:
        print('Running sequentially')
        n_lines = 1
        n_lines_end=number_of_rows_to_load
        for file_in in scored_motifs_overlapping_tracks_files:
            selected_rows_df= pd.read_csv(file_in, nrows=number_of_rows_to_load, sep='\t')
            selected_rows_df_order = selected_rows_df.reindex(columns=col_list[1::])
            #to list of tuples
            records = selected_rows_df_order.to_records(index=False)
            selected_rows = list(records)
        while selected_rows:
            t_to_insert = time.time()
            insert_into_tissues(selected_rows, tissue_cell_assays, tissue_cell_allassays, assay_names,
                                   cols_to_write_to, cols_to_write_to_allassays, thread_num, feature_weights_dict,
                                   db_name, db_user_name, db_host_name, tissues_fscores_table)
            
            print(('t_to_insert: ', time.time()-t_to_insert))
            num_rows -=len(selected_rows)
            print(('num_rows remaining: ', num_rows))
            t_for_fetch = time.time()
            selected_rows = []
            selected_rows_df= pd.read_csv(file_in,skiprows=list(range(n_lines,n_lines_end)), nrows=number_of_rows_to_load, sep='\t', dtype=str)
            selected_rows_df_order = selected_rows_df.reindex(columns=col_list[1::])
            #to list of tuples
            records = selected_rows_df_order.to_records(index=False)
            selected_rows = list(records)
            n_lines_end = n_lines_end + number_of_rows_to_load
            #selected_rows = curs_for_selection.fetchmany(size=number_of_rows_to_load)
            thread_num+=1
            print(('t_to_fetch: ', time.time()-t_for_fetch))
            if num_rows<=0:
                print(("All rows are processed and inserted from {} into tissue tables".format(table_from)))
                break
    
    #if num_rows==0:
    #    print("All rows are processed and inserted from {} into tissue tables".format(table_from))
    #else:
    #    print("Warning!!! There are {} remaining rows to be inserted".format(num_rows))
    print('All rows are processed and inserted into tissues tables')
    #os.remove(tissue_cell_assays_file)
    #os.remove(tissue_cell_allassays_file)
    #os.remove(feature_weights_dict_file)
    #curs_for_selection.close()
    #conn.close()
    return


def populate_tissue_values_from_scored_files_per_tissue(tissue_cell_assays, tissue_cell_allassays, assay_names, col_list, 
                           motif_cols_names,
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
    thread_num = 0
#     conn = DBUtilities.open_connection(db_name, db_user_name, db_host_name)
#     curs_for_count = conn.cursor(name = "countcurs", cursor_factory=DictCursor)
#     thread_num = 0
#     curs_for_count.execute('select count(posrange) from {}'.format(table_from))
#     num_rows = int(curs_for_count.fetchone()[0])#85459976
#     curs_for_count.close()
#     print('total number of rows to be inserted: ', num_rows)
#     curs_for_selection = conn.cursor(name = "selectioncurs", cursor_factory=DictCursor)
#     curs_for_selection.itersize = number_of_rows_to_load
#     t_for_select =  time.time()
#     curs_for_selection.execute('select {} from {}'.format(','.join(col_list), table_from))
    #print('t_to_select: ', time.time() - t_for_select)
    #t_for_fetch = time.time()
    #selected_rows = curs_for_selection.fetchmany(size=number_of_rows_to_load)
    #print('t_to_fetch: ', time.time()-t_for_fetch)
    #print("Selected {} rows for insertion.".format(str(number_of_rows_to_load)))
    list_tissues = list(sorted(tissue_cell_allassays.keys()))
    num_tissues = len(list_tissues)
    t_jobset = time.time()
    if run_in_parallel_param and len(scored_motifs_overlapping_tracks_files)>1:
        print('Running in parallel')
       # n_lines = 1
        #n_lines_end=number_of_rows_to_load
        num_cores = number_processes_to_run_in_parallel
        
        for file_in in scored_motifs_overlapping_tracks_files:
            i = 0
            #n_lines = 1
           # n_lines_end=number_of_rows_to_load
            num_cores = number_processes_to_run_in_parallel
            p = Pool(number_processes_to_run_in_parallel)
            #selected_rows_df= pd.read_csv(file_in, nrows=number_of_rows_to_load, sep='\t', dtype=str)
           # selected_rows_df_order = selected_rows_df.reindex(columns=col_list[1::])
            #to list of tuples
            #records = selected_rows_df_order.to_records(index=False)
            #selected_rows = list(records)
            

            
            
            while i<num_cores:
                
                tissue = list_tissues[i]
                print(tissue)
                scored_motifs_overlapping_tracks_files_tissue = file_in + '_' + tissue
                #cell_assay for tissue
                tissue_cell_assays_v=[]
                tissue_cell_assays_v.extend(motif_cols_names)
                tissue_cell_assays_v = list(tissue_cell_assays[tissue].values())
                #find indeces of tissue-related columns
                tissue_cell_assays_id = [i for i in range(len(col_list)) if col_list[i] in tissue_cell_assays_v]
                #select columns for selected tissue
                comm_cut_tissue = ('cut -f{} {} > {}'.format(','.join([str(id) for id in tissue_cell_assays_id]), file_in, scored_motifs_overlapping_tracks_files_tissue))
                print(comm_cut_tissue)
                os.system(comm_cut_tissue)   
                
                
                p.apply_async(insert_into_tissues_per_tissue, args=(scored_motifs_overlapping_tracks_files_tissue, tissue_cell_assays, tissue_cell_allassays, assay_names,
                                       cols_to_write_to, cols_to_write_to_allassays, thread_num, feature_weights_dict,
                                       db_name, db_user_name, db_host_name, tissues_fscores_table))
                #num_rows -=len(selected_rows)
                print(file_in)
                #print('num_rows: ', n_lines_end)
                
                #selected_rows = []
                #selected_rows_df= pd.read_csv(file_in,skiprows=range(n_lines,n_lines_end), nrows=number_of_rows_to_load, sep='\t', dtype=str)
                #selected_rows_df_order = selected_rows_df.reindex(columns=col_list[1::])
                #to list of tuples
                #records = selected_rows_df_order.to_records(index=False)
                #selected_rows = list(records)
                #n_lines_end = n_lines_end + number_of_rows_to_load
                #selected_rows = curs_for_selection.fetchmany(size=number_of_rows_to_load)
                num_tissues -=1
                if num_tissues<=0:
                    p.close()
                    p.join()
                    break
                if i==num_cores-1:
                    p.close()
                    p.join()
                    print(('t_jobset: ', time.time()-t_jobset))
                    t_jobset = time.time()
                    p = Pool(number_processes_to_run_in_parallel)
                    i=0
                i+=1
                thread_num+=1
    else:
        print('Running sequentially')
        n_lines = 1
        n_lines_end=number_of_rows_to_load
        for file_in in scored_motifs_overlapping_tracks_files:
            selected_rows_df= pd.read_csv(file_in, nrows=number_of_rows_to_load, sep='\t')
            selected_rows_df_order = selected_rows_df.reindex(columns=col_list[1::])
            #to list of tuples
            records = selected_rows_df_order.to_records(index=False)
            selected_rows = list(records)
        while selected_rows:
            t_to_insert = time.time()
            insert_into_tissues(selected_rows, tissue_cell_assays, tissue_cell_allassays, assay_names,
                                   cols_to_write_to, cols_to_write_to_allassays, thread_num, feature_weights_dict,
                                   db_name, db_user_name, db_host_name, tissues_fscores_table)
            
            print(('t_to_insert: ', time.time()-t_to_insert))
            num_rows -=len(selected_rows)
            print(('num_rows remaining: ', num_rows))
            t_for_fetch = time.time()
            selected_rows = []
            selected_rows_df= pd.read_csv(file_in,skiprows=list(range(n_lines,n_lines_end)), nrows=number_of_rows_to_load, sep='\t', dtype=str)
            selected_rows_df_order = selected_rows_df.reindex(columns=col_list[1::])
            #to list of tuples
            records = selected_rows_df_order.to_records(index=False)
            selected_rows = list(records)
            n_lines_end = n_lines_end + number_of_rows_to_load
            #selected_rows = curs_for_selection.fetchmany(size=number_of_rows_to_load)
            thread_num+=1
            print(('t_to_fetch: ', time.time()-t_for_fetch))
            if num_rows<=0:
                print(("All rows are processed and inserted from {} into tissue tables".format(table_from)))
                break
    
    if num_rows==0:
        print(("All rows are processed and inserted from {} into tissue tables".format(table_from)))
    else:
        print(("Warning!!! There are {} remaining rows to be inserted".format(num_rows)))
    #curs_for_selection.close()
    #conn.close()
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
                    print("no weights is used for sl[0] because the assigned values is not a number")
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
                    generate_tissue_from_db=False
                    ):
    col_list, tissue_cell_assays, tissue_cell_allassays = get_tissue_cell_mappings(cell_assays, 
                                                                                       assay_names, 
                                                                                       tissue_cell_mappings_file,  
                                                                                       motif_cols = motif_cols_names)
    print(col_list)
    tissue_cols = db_setup_tissues(tissue_cell_allassays, 
                                   assay_cells_datatypes, 
                                   motif_cols = ['mid INTEGER', 'fscore NUMERIC'], 
                                   db_name = db_name, db_user_name=db_user_name, db_host_name=db_host_name,
                                   tissues_fscores_table=tissues_fscores_table,
                                   motif_cols_for_tissues_fscores_table = ['mid INTEGER'])
    col_list.append('mid')
    feature_weights_dict = get_weights_per_feature(annotation_weights_inputfile=annotation_weights_inputfile,skip_negative_weights=skip_negative_weights)
    print('Inserting data into tissues tables')
    if generate_tissue_from_db:
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
    else:
        populate_tissue_values_from_scored_files(tissue_cell_assays, 
                               tissue_cell_allassays, 
                               assay_names, col_list, 

                               run_in_parallel_param=run_in_parallel_param, 
                               number_processes_to_run_in_parallel=number_processes_to_run_in_parallel, 
                               scored_motifs_overlapping_tracks_files=scored_motifs_overlapping_tracks_files, 
                               feature_weights_dict=feature_weights_dict,
                               db_name = db_name, db_user_name=db_user_name, db_host_name=db_host_name,
                               tissues_fscores_table=tissues_fscores_table,
                               number_of_rows_to_load = number_of_rows_to_load,
                               )
        
        
    
    print("Creating index on tissues tables")
    for tissue_table in sorted(tissue_cols.keys()):
        DBUtilities.create_index(db_name, db_user_name, db_host_name, 
                                 cell_table=tissue_table, index_name='index'+tissue_table+'mid', index_method = 'btree', index_cols = 'mid')
    DBUtilities.create_index(db_name, db_user_name, db_host_name, 
                                 cell_table=tissues_fscores_table, index_name='index'+tissues_fscores_table+'mid', index_method = 'btree', index_cols = 'mid')
