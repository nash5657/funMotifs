import psycopg2

# Your database connection settings
database_name = 'funmotifsdb'
user_name = 'naser'
host_address = '127.0.0.1'

# List of tables to query
tables = ['blood', 'brain', 'breast', 'cervix', 'colon', 'esophagus', 'kidney', 'liver', 'lung', 'myeloid', 'pancreas', 'prostate', 'skin', 'stomach', 'uterus']

output_file = 'numeric_values_all_tables.txt'
try:
    # Establishing the connection
    conn = psycopg2.connect(
        database=database_name,
        user=user_name,
        host=host_address
    )

    # Creating a cursor object
    cursor = conn.cursor()
    # List of columns to check
    columns_to_check = ['mid', 'fscore', 'contactingdomain', 'dnase__seq', 'fantom', 'numothertfbinding', 'tfbinding', 'tfexpr', 'footprints']
    # Writing results to a file
    with open(output_file, 'w') as file:
        for table in tables:
            file.write(f"Table: {table}\n")
            for column in columns_to_check:
                # SQL query for each column in each table
                query = f"SELECT EXISTS(SELECT 1 FROM {table} WHERE {column} != 0 LIMIT 1);"
                # Executing the query
                cursor.execute(query)

                # Fetching the results
                results = cursor.fetchall()

                # Writing column name and its values
                file.write(f"    Column: {column}\n")
                for row in results:
                    file.write(f"        {row[0]}\n")
                file.write("\n")  # Adding a newline for separation
            file.write("\n")  # Separating different tables
except Exception as e:
    print(f"An error occurred: {e}")
finally:
    # Closing the cursor and connection
    if cursor:
        cursor.close()
    if conn:
        conn.close()
