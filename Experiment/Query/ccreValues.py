import psycopg2

# Your database connection settings
database_name = 'funmotifsdb'
user_name = 'naser'
host_address = '127.0.0.1'

# List of tables to query
tables = ['blood', 'brain', 'breast', 'cervix', 'colon', 'esophagus', 'kidney', 'liver', 'lung', 'myeloid', 'pancreas', 'prostate', 'skin', 'stomach', 'uterus']

output_file = 'ccre_values_all_tables.txt'

try:
    # Establishing the connection
    conn = psycopg2.connect(
        database=database_name,
        user=user_name,
        host=host_address
    )

    # Creating a cursor object
    cursor = conn.cursor()

    # Open the output file
    with open(output_file, 'w') as file:
        for table in tables:
            # SQL query
            query = f"SELECT DISTINCT ccre FROM {table} WHERE ccre IS NOT NULL;"

            # Executing the query
            cursor.execute(query)

            # Fetching all results
            results = cursor.fetchall()

            # Writing results to the file
            for row in results:
                file.write(f"{table}, {row[0]}\n")

except Exception as e:
    print(f"An error occurred: {e}")
finally:
    # Closing the cursor and connection
    if cursor:
        cursor.close()
    if conn:
        conn.close()
