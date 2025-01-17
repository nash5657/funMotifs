'''
Created on 6 Oct 2017

@author: husensofteng
@contributor: markmelzer
@contributer: Naser Shabani
'''
import sys
import os
import psycopg2
from psycopg2.extensions import ISOLATION_LEVEL_AUTOCOMMIT


def start_psql_server(db_user_name, db_host_name, db="../database", logfile="../logfile"):
    # start the PostgreSQL server
    try:
        # test if server is already running
        con_postgres = psycopg2.connect(dbname='postgres', user=db_user_name, host=db_host_name)
        curs = con_postgres.cursor()
        curs.close()
        con_postgres.close()

    except:
        if sys.platform == "linux" or sys.platform == "linux2":
            try:
                if os.path.exists(db):
                    os.system("""pg_ctl restart -D """ + db + """ -l """ + logfile)
                else:
                    raise Exception("Path does not yet exist")
            except:
                try:
                    if os.path.exists(db):
                        os.system("""pg_ctl initdb -D """ + db)
                        os.system("""pg_ctl -D """ + db + """ -l logfile start""")
                    else:
                        os.system("""mkdir """ + db)
                        os.system("""pg_ctl initdb -D """ + db)
                        os.system("""pg_ctl -D """ + db + """ -l logfile start""")
                except:
                    raise Exception("PostgreSQL server could not be started. Please check your permissions or start "
                                    "database  manually, then run the program again.")
        else:
            raise Exception("First start the PostgreSQL server, then run the program again.")
    return


def create_db(db_name, db_user_name, db_host_name):
    conn = ""
    curs = ""
    try:
        conn = psycopg2.connect("dbname={} user={} host={}".format(db_name, db_user_name, db_host_name))
        curs = conn.cursor()
        curs.execute("SELECT exists(SELECT 1 from pg_catalog.pg_database where datname = %s)", (db_name,))
        if curs.fetchone()[0]:
            print(("Successfully connected to DB: ", db_name))

    except psycopg2.DatabaseError as e:
        print(("Error", e))
        print(("Creating DB: ", db_name))
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
            print(("Successfully created and connected to DB: ", db_name))
            conn.close()
            return True
        else:
            return False



def open_connection(db_name, db_user_name='user', db_host_name='localhost'):
    conn = psycopg2.connect("dbname={} user={} host={}".format(db_name, db_user_name, db_host_name))
    # conn.row_factory = psycopg2.Row
    return conn


def close_connection(conn):
    conn.close()


def get_col_names_from_table(table_name, conn):
    curs = conn.cursor()
    curs.execute("select * FROM {} limit 1".format(table_name))
    return [desc[0] for desc in curs.description]

def create_index(cursor, cell_table, index_name='indexposrange', index_method='gist', index_cols='posrange'):
    """
    Function to create an index in the database. Uses an existing database cursor.
    """
    cursor.execute(f"DROP INDEX IF EXISTS {index_name}")
    create_index_stmt = f"CREATE INDEX IF NOT EXISTS {index_name} ON {cell_table} USING {index_method} ({index_cols})"
    # print(create_index_stmt)
    cursor.execute(create_index_stmt)
    cursor.connection.commit()

    
# def create_index(db_name, db_user_name, db_host_name,
#                  cell_table, index_name='indexposrange', index_method='gist', index_cols='posrange'):
#     conn = open_connection(db_name, db_user_name, db_host_name)

#     curs = conn.cursor()
#     curs.execute("DROP INDEX IF EXISTS {}".format(index_name))
#     creat_index_stmt = "CREATE INDEX IF NOT EXISTS {} ON {} using {} ({})".format(index_name, cell_table, index_method,
#                                                                                   index_cols)
#     # print creat_index_stmt
#     curs.execute(creat_index_stmt)
#     conn.commit()
#     close_connection(conn)
#     return

def table_contains_data(cursor, table_name):
    """
    Function to check if a table contains data. Uses an existing database cursor.
    """
    try:
        cursor.execute(f"SELECT chr FROM {table_name} LIMIT 1")
        if cursor.fetchone() is not None:
            print(f"{table_name} contains data")
            return True
        else:
            return False
    except psycopg2.ProgrammingError as e:
        print(f"Error: {e}")
        return False
    
# def table_contains_data(db_name, db_user_name, db_host_name, table_name):
#     conn = open_connection(db_name, db_user_name, db_host_name)
#     curs = conn.cursor()
#     try:
#         curs.execute('select chr from {} limit 1'.format(table_name))
#         if curs.fetchone() is not None:
#             print(('{} contains data'.format(table_name)))
#             return True
#         else:
#             return False
#     except psycopg2.ProgrammingError:
#         return False
#     close_connection(conn)

def create_table_stmt_parallel(cursor, tissue, tissuecols, tissuemotifsimputed):
    """
    Function to create a table in the database. Uses an existing database cursor.
    """
    
    create_table_query = f"CREATE TABLE IF NOT EXISTS {tissue} AS (SELECT {tissuecols} FROM {tissuemotifsimputed})"
    cursor.execute(create_table_query)
    cursor.connection.commit()


# def create_table_stmt_parallel(db_name, db_user_name, db_host_name,
#                                tissue, tissuecols, tissuemotifsimputed):
#     conn = open_connection(db_name, db_user_name, db_host_name)
#     curs = conn.cursor()
#     # print "Loading for", tissue
#     # print tissuecols
#     # print "create table if not exists {0} as (select {1} from {2})".format(tissue, tissuecols, tissuemotifsimputed)
#     curs.execute(
#         "create table if not exists {0} as (select {1} from {2})".format(tissue, tissuecols, tissuemotifsimputed))
#     conn.commit()
#     # print "Created", tissue
#     curs.close()
#     conn.close()

def get_number_of_motifs(table: str, cursor):
    """
    Function that returns the total number of motifs saved in the table.
    Now takes a database cursor as an argument.
    """
    cursor.execute(f"SELECT count(*) FROM {table}")
    num = cursor.fetchone()[0]
    return num

# def get_number_of_motifs(table: str, db_name: str, db_user_name: str, db_host_name: str):
#     """
#     Function that returns the total number of motifs saved in the table
#     """
#     #conn = psycopg2.connect(database=db_name, user=db_user_name)
#     conn = open_connection(db_name, db_user_name, db_host_name)
#     curs = conn.cursor()
#     curs.execute(f"""SELECT count(*) from {table}""")
#     num = curs.fetchone()[0]
#     conn.commit()
#     curs.close()
#     return num

def add_column_to_tissue_table(table: str, cursor, col_name: str, col_type: str):
    """
    Function adds an additional column name to the annotated motif data of a tissue.
    Now takes a database cursor as an argument.
    """
    cursor.execute(f"ALTER TABLE {table} ADD COLUMN IF NOT EXISTS {col_name} {col_type}")
    cursor.connection.commit()

# def add_column_to_tissue_table(table: str, cursor, db_name: str, db_user_name: str, db_host_name: str, col_name: str, col_type: str):
#     """
#     Function adds an additional column name to the annotated motif data of a tissue
#     """
    # #conn = psycopg2.connect(database=db_name, user=db_user_name)
    # conn = open_connection(db_name, db_user_name, db_host_name)
    # cur = conn.cursor()
    # cur.execute(f"ALTER TABLE {table} ADD COLUMN IF NOT EXISTS {col_name} {col_type}")
    # conn.commit()
    # cur.close()
    # return

def update_db_value(value, mid, column, table, cursor):
    """
    Function to update a value in a Postgres database.
    Now takes a database cursor as an argument.
    """
    sql = f"UPDATE {table} SET {column} = {value} WHERE mid = {mid}"
    cursor.execute(sql)
    cursor.connection.commit()

# def update_db_value(value, mid, column, table, db_name, db_user_name):
#     """
#     Function to update a value in a Postgres database
#     """
#     sql = f"UPDATE {table} SET {column} = {value} WHERE mid = {mid}"
#     conn = psycopg2.connect(database=db_name, db_user_name=db_user_name)
#     cur = conn.cursor()
#     cur.execute(sql)
#     conn.commit()
#     cur.close()
#     return

