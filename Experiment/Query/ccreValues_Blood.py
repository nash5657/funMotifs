
import psycopg2

# Your database connection settings
database_name = 'funmotifsdb'
user_name = 'naser'
host_address = '127.0.0.1'

output_file = 'ccre_values.txt'

try:
    # Establishing the connection
    conn = psycopg2.connect(
        database=database_name,
        user=user_name,
        host=host_address
    )

    # Creating a cursor object
    cursor = conn.cursor()

    # SQL query
    query = "SELECT DISTINCT ccre FROM blood WHERE ccre IS NOT NULL AND ccre <> 'NO';"

    # Executing the query
    cursor.execute(query)

    # Fetching the results
    results = cursor.fetchall()

    # Writing results to a file
    with open(output_file, 'w') as file:
        for row in results:
            file.write(f"{row[0]}
")

    # Closing the cursor and connection
    cursor.close()
    conn.close()

    print(f"Results saved to {output_file}")

except Exception as e:
    print(f"An error occurred: {e}")
