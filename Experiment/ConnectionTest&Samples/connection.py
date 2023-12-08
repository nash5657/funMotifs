import psycopg2

# Replace with your database information
database_name = 'funmotifsdb'
user_name = 'naser'
host_address = '127.0.0.1'
# password = 'your_password'
# port = '5432'

output_file = 'database_structure_and_data.txt'

try:
    # Establishing the connection
    conn = psycopg2.connect(
        database=database_name,
        user=user_name,
        host=host_address
        # password=password,
        # port=port
    )

    # Creating a cursor object
    cursor = conn.cursor()

    # Opening the file to write
    with open(output_file, 'w') as file:
        # Query to get all table names
        cursor.execute("SELECT table_name FROM information_schema.tables WHERE table_schema = 'public'")
        tables = cursor.fetchall()

        for table in tables:
            table_name = table[0]

            # Query for table structure
            cursor.execute(f"SELECT column_name, data_type FROM information_schema.columns WHERE table_name = '{table_name}'")
            columns = cursor.fetchall()
            file.write(f"\nStructure of table '{table_name}':\n")
            for column in columns:
                file.write(f"{column[0]} ({column[1]})\n")

            # # Query to get row count for the table
            # cursor.execute(f"SELECT COUNT(*) FROM {table_name}")
            # row_count = cursor.fetchone()[0]
            # file.write(f"Number of rows in '{table_name}': {row_count}\n")

            # Query for sample data
            cursor.execute(f"SELECT * FROM {table_name} LIMIT 5")
            sample_data = cursor.fetchall()
            file.write(f"Sample data from '{table_name}':\n")
            for row in sample_data:
                file.write(f"{row}\n")

except Exception as e:
    print("An error occurred:", e)
finally:
    # Closing the cursor and connection
    if cursor:
        cursor.close()
    if conn:
        conn.close()

print(f"Results have been saved to {output_file}")
