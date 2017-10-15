'''
Created on 6 Oct 2017

@author: husensofteng
'''
from multiprocessing import Pool
import psycopg2
from psycopg2.extensions import ISOLATION_LEVEL_AUTOCOMMIT


def create_db(db_name, db_user_name, db_host_name):
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
            conn.close()
            return True
        else:
            return False

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


def create_index(db_name, db_user_name, db_host_name, 
                 cell_table, index_name='indexposrange', index_method = 'gist', index_cols = 'posrange'):
    conn = open_connection(db_name, db_user_name, db_host_name)
        
    curs = conn.cursor()
    curs.execute("DROP INDEX IF EXISTS {}".format(index_name))
    creat_index_stmt= "CREATE INDEX IF NOT EXISTS {} ON {} using {} ({})".format(index_name, cell_table, index_method, index_cols)
    print creat_index_stmt
    curs.execute(creat_index_stmt)
    conn.commit()
    close_connection(conn)
    return

def table_contains_data(db_name, db_user_name, db_host_name, table_name):
    conn = open_connection(db_name, db_user_name, db_host_name)
    curs = conn.cursor()
    try:
        curs.execute('select chr from {} limit 1'.format(table_name))
        if curs.fetchone() is not None:
            print '{} contains data'.format(table_name)
            return True
        else:
            return False
    except psycopg2.ProgrammingError:
        return False
    close_connection(conn)

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

def create_motifs_table(db_name, db_user_name, db_host_name, motifs_table, motif_cols, new_table_name):
    conn = open_connection(db_name, db_user_name, db_host_name)
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
    conn = open_connection(db_name, db_user_name, db_host_name)
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
    print 'All tables are created'
    return

