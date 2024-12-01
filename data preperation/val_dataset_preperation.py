import pandas as pd

# Load the validation dataset
file_path = pd.read_excel("/Users/oliver/PycharmProjects/pythonProject/internship_data/10.xlsx", skiprows=1)

# Select columns
val_data = file_path[['Spacer_ID', 'Spacer Sequence', 'Normalized sgRNA Activity']]

# Rename columns
val_data.columns = ['sgRNAID', 'nucleotide sequence', 'score']

# Save validation data
val_file_path = '/Users/oliver/PycharmProjects/pythonProject/internship_data/val_data.csv'
val_data.to_csv(val_file_path, index=False)
