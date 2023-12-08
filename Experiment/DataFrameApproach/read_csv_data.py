import pandas as pd

# Path to the CSV file - replace with the actual path of your CSV file
csv_file_path = 'blood.csv'  # Replace with your actual file path

# Read the CSV file into a DataFrame
df = pd.read_csv(csv_file_path)

# Display the shape of the DataFrame
print("DataFrame shape:", df.shape)

# Display the first few rows of the DataFrame as a sample
print(df.head())
