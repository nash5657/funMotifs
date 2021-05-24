'''
Created on 6 Oct 2017

@author: husensofteng
'''
import sys
from itertools import islice
from multiprocessing import Pool
import psycopg2
import DBUtilities

def create_cell_table(db_name, db_user_name, db_host_name,
                      cells_assays_dict, 
                         assay_cells_datatypes, 
                         cell_table, 
                         motif_col_names, 
                         motif_cols,
                         ):
    
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
            
            if cell=='22Rv1' or cell=='8988T':
                cell='a'+cell
                
            if cell=="Ammon's horn":
                cell="Ammons horn"
            if cell=="Peyer's patch":
                cell="Peyers patch"

                

            feature ='_'.join(((cell + "___" + assay).replace('(','').replace(')','')
                                         .replace('-','__').replace('.','')).split())
            field_names.append(feature+ " " + data_type)
            
            col_names.append(feature)
    print(field_names)
    #curs.execute("DROP TABLE IF EXISTS {}".format(cell_table))
    conn = DBUtilities.open_connection(db_name, db_user_name, db_host_name)
    curs = conn.cursor()
    create_table_stmt = "CREATE TABLE IF NOT EXISTS {} ({});".format(cell_table, ', '.join(field_names))
    curs.execute(create_table_stmt)
    conn.commit()
    DBUtilities.close_connection(conn)
    return col_names


def insert_from_file(field_names, i_file, n, db_name, db_user_name, 
                     db_host_name, cell_table, header,
                     thread_num=0):
    
    n_processed = 0
    conn = DBUtilities.open_connection(db_name, db_user_name, db_host_name)
    curs = conn.cursor()
    with open(i_file, 'r') as ifile:
        print(i_file)
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
                    print("Thread {} for ({}) has processed: {}".format(thread_num, i_file, n_processed))
                except:
                    print("Thread {} coundn't insert to DB ({})".format(thread_num, i_file))
            except (ValueError, IndexError):
                print("Empty file" + i_file + "in thread" + str(thread_num))
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
    log_file = open("log_file.txt", 'w')
    if run_in_parallel_param and len(scored_motifs_overlapping_tracks_files)>1:
        p = Pool(number_processes_to_run_in_parallel)
    thread_num=0
    for i_file in scored_motifs_overlapping_tracks_files:#[f for f in glob.glob('{}/*{}*'.format(dir_to_import, keyword_to_check))]
        if run_in_parallel_param and len(scored_motifs_overlapping_tracks_files)>1:
            thread_num+=1
            p.apply_async(insert_from_file, args=[field_names, i_file, 100000, db_name, db_user_name, db_host_name, cell_table, header, thread_num])
        else:
            insert_from_file(field_names, i_file, 100000, db_name, db_user_name, db_host_name, cell_table, header)
        log_file.write(i_file+'\n')
    if run_in_parallel_param and len(scored_motifs_overlapping_tracks_files)>1:
        p.close()
        p.join()
    
    print("Data insertion into {} is done".format(cell_table))
    return

def generate_cell_table(db_name,
                    cell_table,
                    db_user_name,
                    db_host_name,
                    cells_assays_dict,
                    assay_cells_datatypes,
                    run_in_parallel_param,
                    number_processes_to_run_in_parallel,
                    header,
                    scored_motifs_overlapping_tracks_files,
                    motif_cols,
                    motif_cols_names,
                    cell_index_name, cell_index_method, cell_index_cols
        ):
    
    if not DBUtilities.table_contains_data(db_name, db_user_name, db_host_name, 
                                           cell_table): #that is to avoid writing over an already existing table content
        
        field_names = create_cell_table(db_name, db_user_name, db_host_name,
            cells_assays_dict, assay_cells_datatypes, cell_table, motif_cols_names[1:], motif_cols)
        
        print("Inserting data into: ", cell_table)
        insert_into_db(field_names, db_name = db_name, 
                       db_user_name=db_user_name, 
                       db_host_name=db_host_name, 
                       cell_table=cell_table, 
                       scored_motifs_overlapping_tracks_files=scored_motifs_overlapping_tracks_files, 
                       header=header,
                       run_in_parallel_param=run_in_parallel_param,
                       number_processes_to_run_in_parallel=number_processes_to_run_in_parallel)
        # dir_to_import=params['motifs_overlapping_tracks_output_dir'], keyword_to_check="_scored.bed10", header=header)
        print("Creating index on: ", cell_table)
        DBUtilities.create_index(db_name, db_user_name, db_host_name, 
                                 cell_table, index_name=cell_index_name, index_method = cell_index_method, index_cols = cell_index_cols)
        
