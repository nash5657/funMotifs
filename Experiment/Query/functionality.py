import psycopg2

# Replace with your database information
database_name = 'funmotifsdb'
user_name = 'naser'
host_address = '127.0.0.1'
# password = 'your_password'
# port = '5432'

output_file = 'functional_motifs.txt'

try:
    # Establishing the connection
    conn = psycopg2.connect(
        database=database_name,
        user=user_name,
        host=host_address,
        # password=password,
        # port=port
    )

    # Creating a cursor object
    cursor = conn.cursor()

    # SQL Query to select motifs with non-null functionality
    sql_query = "SELECT * FROM motifs WHERE functionality IS NOT NULL;"

    # Execute the query
    cursor.execute(sql_query)

    # Fetch all the rows
    rows = cursor.fetchall()

    # Write to a text file
    with open(output_file, 'w') as file:
        for row in rows:
            file.write(str(row) + '\\n')

    # Closing the cursor and connection
    cursor.close()
    conn.close()

except Exception as e:
    print(f"Database connection failed due to {e}")

