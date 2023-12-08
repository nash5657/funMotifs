import psycopg2

# Database connection details
database_name = 'funmotifsdb'
user_name = 'naser'
host_address = '127.0.0.1'
# Uncomment and set your password and port if necessary
# password = 'your_password'
# port = '5432'

# Establishing the connection
conn = psycopg2.connect(
    database=database_name,
    user=user_name,
    host=host_address,
    # password=password,
    # port=port
)

# Define your SQL query - replace with your actual query
sql_query = "COPY (SELECT * FROM blood) TO STDOUT WITH CSV HEADER;"  # Replace 'your_table_name' with the actual table name

# Path to the output CSV file
output_csv_file = 'blood.csv'  # Replace with your desired file path

try:
    with conn.cursor() as cursor:
        with open(output_csv_file, 'w') as f_output:
            cursor.copy_expert(sql_query, f_output)
finally:
    # Close the database connection
    conn.close()

print(f"Query results saved to {output_csv_file}")
