import pandas as pd

# load data from each sheet in Data S1.xlsx
file_path_s1 = "/Users/oliver/PycharmProjects/pythonProject/internship_data/Data S1.xlsx"
promoter = pd.read_excel(file_path_s1, sheet_name='promoter')
rbs = pd.read_excel(file_path_s1, sheet_name='RBS')
crispri = pd.read_excel(file_path_s1, sheet_name='CRISPRi (gene-targeting)')

# combine data frames into one
sgRNA_sequence = pd.concat([promoter, rbs, crispri], ignore_index=True)
print(sgRNA_sequence.head())

# save sgRNA sequence
sgRNA_file_path = '/Users/oliver/PycharmProjects/pythonProject/internship_data/sgRNA_sequence.csv'
sgRNA_sequence.to_csv(sgRNA_file_path, index=False)

# load Cas9 data from Data S6.xlsx
file_path_s6 = "/Users/oliver/PycharmProjects/pythonProject/internship_data/Data S6.xlsx"
data_s6 = pd.read_excel(file_path_s6, sheet_name='Cas9')
print(data_s6.head())

# merge two files
merged_data = pd.merge(sgRNA_sequence, data_s6, left_on='sgRNAID', right_on='sgRNA', how='inner')
print(merged_data.head())

# save merged data
merged_file_path = '/Users/oliver/PycharmProjects/pythonProject/internship_data/merged_data.csv'
merged_data.to_csv(merged_file_path, index=False)

# Rename Log2_normalized_change to score and round it to second decimal
merged_data.rename(columns={'Log2_normalized_change': 'score'}, inplace=True)
merged_data['score'] = merged_data['score'].round(2)

# drop columns sgRNA & gene (promoter)
filtered_data = merged_data.drop(columns=['sgRNA', 'gene (promoter)'])

# drop Quality = Bad
ecoli_data = filtered_data[filtered_data['Quality'] != 'Bad']
print(ecoli_data.head())
row_count = ecoli_data.shape[0]
print(f"ecoli_data has {row_count} rows of data")

# save ecoli data
ecoli_file_path = '/Users/oliver/PycharmProjects/pythonProject/internship_data/ecoli_data.csv'
ecoli_data.to_csv(ecoli_file_path, index=False)
