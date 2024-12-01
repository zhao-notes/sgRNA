import pandas as pd

onehot_dict = {
    'A': '1000',
    'C': '0100',
    'T': '0010',
    'G': '0001'
}

# Generate the feature names from the dictionary
feature_names = [f'{key}_{i + 1}' for key in onehot_dict.keys() for i in range(20)]

# Input and output file paths
input_path = '/Users/oliver/PycharmProjects/pythonProject/internship_data/val_data.csv'
output_path = '/Users/oliver/PycharmProjects/pythonProject/internship_data/val_dep1.txt'

# Open the input file and read it into a DataFrame
df = pd.read_csv(input_path)

# Initialize a list to store the rows for the output DataFrame
output_rows = []

# Loop over nucleotide sequences
for idx, row in df.iterrows():
    sgRNAID = row['sgRNAID']
    seq = row['nucleotide sequence']

    # Initialize a dictionary to store the feature values
    feature_dict = {feature: 0 for feature in feature_names}

    # Set sgRNAID as the first key-value pair in the dictionary
    feature_dict = {'sgRNAID': sgRNAID, **feature_dict}

    # Compute position-dependent features and update the dictionary
    for i in range(len(seq)):
        monomer = seq[i]
        feature_name = f'{monomer}_{i + 1}'
        if feature_name in feature_dict:
            feature_dict[feature_name] = 1

    # Append the dictionary to the list of rows
    output_rows.append(feature_dict)

# Convert the list of rows to a DataFrame
output_df = pd.DataFrame(output_rows)

# Save the output DataFrame to a text file
output_df.to_csv(output_path, sep=',', index=False)

print(f'Data saved to {output_path}')

