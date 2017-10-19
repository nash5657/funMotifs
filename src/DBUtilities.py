'''
Created on 6 Oct 2017

@author: husensofteng
'''
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
    conn.commit()
    print "Created", tissue
    curs.close()
    conn.close()
