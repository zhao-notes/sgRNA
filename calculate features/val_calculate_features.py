import pandas as pd
import numpy as np
import os, sys
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import GC
from ViennaRNA import RNA

# Load validation data
file_path = '/Users/oliver/PycharmProjects/pythonProject/internship_data/val_data.csv'
df = pd.read_csv(file_path, low_memory=False)
df.info()

# Extract sequences for calculations
sequences = df['nucleotide sequence']

# Calculate Tm and GC content
df['Tm'] = sequences.apply(mt.Tm_Wallace)
df['GC_content'] = sequences.apply(GC)


# Function to calculate Minimum Free Energy (MFE) using ViennaRNA
def calculate_mfe(sequence):
    (mfe_structure, mfe) = RNA.fold(sequence)
    return mfe


df['MFE'] = sequences.apply(calculate_mfe).round(2)


# Function to compute position-independent features
def compute_pos_indep_features(seq):
    # Single nucleotide counts
    single_counts = {nuc: seq.count(nuc) for nuc in 'ACGT'}
    # Paired nucleotide counts
    pair_counts = {f"{n1}{n2}": seq.count(f"{n1}{n2}") for n1 in 'ACGT' for n2 in 'ACGT'}
    return single_counts, pair_counts

# Compute position-independent features and add to dataframe
single_features = []
pair_features = []

for seq in sequences:
    single_counts, pair_counts = compute_pos_indep_features(seq)
    single_features.append(single_counts)
    pair_features.append(pair_counts)

single_df = pd.DataFrame(single_features)
pair_df = pd.DataFrame(pair_features)

# Add sgRNAID
single_df['sgRNAID'] = df['sgRNAID']
pair_df['sgRNAID'] = df['sgRNAID']

# Reorder columns to keep 'sgRNAID' as the first column
single_df = single_df[['sgRNAID'] + [col for col in single_df.columns if col != 'sgRNAID']]
pair_df = pair_df[['sgRNAID'] + [col for col in pair_df.columns if col != 'sgRNAID']]

# Save single and pair features to separate text files
single_output_path = '/Users/oliver/PycharmProjects/pythonProject/internship_data/val_ind1.txt'
pair_output_path = '/Users/oliver/PycharmProjects/pythonProject/internship_data/val_ind2.txt'

single_df.to_csv(single_output_path, sep=',', index=False)
pair_df.to_csv(pair_output_path, sep=',', index=False)

print(f"Single base features saved to {single_output_path}")
print(f"Paired base features saved to {pair_output_path}")



# Prepare the data for MFE.txt
mfe_output_path = '/Users/oliver/PycharmProjects/pythonProject/internship_data/val_MFE.txt'
mfe_df = df[['sgRNAID', 'MFE']]
mfe_df.to_csv(mfe_output_path, sep='\t', index=False)
print(f"MFE data saved to {mfe_output_path}")

# Prepare the data for nuc.count.txt
nuc_output_path = '/Users/oliver/PycharmProjects/pythonProject/internship_data/val_nuc.count.txt'
nuc_df = df[['sgRNAID', 'Tm', 'GC_content']]
nuc_df.to_csv(nuc_output_path, sep='\t', index=False)
print(f"Nucleotide count data saved to {nuc_output_path}")

# Prepare the data for score.txt
score_output_path = '/Users/oliver/PycharmProjects/pythonProject/internship_data/val_score.txt'
score_df = df[['sgRNAID', 'score', 'nucleotide sequence']]  # Ensure 'cutting score' is the correct column name
score_df.to_csv(score_output_path, sep='\t', index=False)
print(f"Score data saved to {score_output_path}")
