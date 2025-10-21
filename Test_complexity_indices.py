import pandas as pd
from testing.diversity_and_complexity.molecule_complexity_indexs import bertz_index, hann_index, wiener_index, sas_index, qed_index

# Read the CSV file
folder_path = "/home/raul-acosta/GitHub/mutation_and_recombination/results/computational_resources/Pairs.csv"
df = pd.read_csv(folder_path)

# Drop rows with NaN values in the 'Mutated_SMILES' column
df = df.dropna(subset=['Mutated_SMILES'])

# Reset the index after dropping rows
df.reset_index(drop=True, inplace=True)

# Function to calculate complexity indices for each SMILES
def calculate_complexity_indices(smiles):
    try:
        bertz_idx = bertz_index(smiles)
        hann_idx = hann_index(smiles)
        wiener_idx = wiener_index(smiles)
        sas_idx = sas_index(smiles)
        qed_idx = qed_index(smiles)
        return bertz_idx, hann_idx, wiener_idx, sas_idx, qed_idx
    except ValueError as e:
        print(f"Error calculating complexity indices for SMILES '{smiles}': {e}")
        return None, None, None, None, None


"""

# Apply the function to calculate complexity indices for each Original_SMILES
df[['Original_Bertz_index', 'Original_Hann_index', 
    'Original_Wiener_index', 'Original_SAS_index', 
    'Original_QED_index']] = pd.DataFrame(df['Original_SMILES'].apply(calculate_complexity_indices).tolist(), index=df.index)

# Apply the function to calculate complexity indices for each Mutated_SMILES
df[['Mutated_Bertz_index', 'Mutated_Hann_index', 
    'Mutated_Wiener_index', 'Mutated_SAS_index', 
    'Mutated_QED_index']] = pd.DataFrame(df['Mutated_SMILES'].apply(calculate_complexity_indices).tolist(), index=df.index)

# Specify the new columns for cleaning
new_columns_to_check = ['Mutated_SMILES', 'Mutated_Bertz_index', 'Mutated_Hann_index', 
                        'Mutated_Wiener_index', 'Mutated_SAS_index', 'Mutated_QED_index']

# Drop rows with NaN or None values in the new specified columns
df = df.dropna(subset=new_columns_to_check)

# Reset the index after dropping rows
df.reset_index(drop=True, inplace=True)

# Save the cleaned DataFrame to a new CSV file
df.to_csv("Pairs_cleaned_complexity_indices.csv", index=False)

"""
df = pd.read_csv("Pairs_cleaned_complexity_indices.csv")

# Calculate the absolute differences for each complexity index
df['Abs_Diff_Bertz_index'] = abs(df['Original_Bertz_index'] - df['Mutated_Bertz_index'])
df['Abs_Diff_Hann_index'] = abs(df['Original_Hann_index'] - df['Mutated_Hann_index'])
df['Abs_Diff_Wiener_index'] = abs(df['Original_Wiener_index'] - df['Mutated_Wiener_index'])
df['Abs_Diff_SAS_index'] = abs(df['Original_SAS_index'] - df['Mutated_SAS_index'])
df['Abs_Diff_QED_index'] = abs(df['Original_QED_index'] - df['Mutated_QED_index'])

# Save the DataFrame with absolute differences to a new CSV file
df.to_csv("Pairs_cleaned_complexity_indices_difference.csv", index=False)

# Calculate the absolute differences for each complexity index
df['Abs_Diff_Bertz_index'] = abs(df['Original_Bertz_index'] - df['Mutated_Bertz_index'])
df['Abs_Diff_Hann_index'] = abs(df['Original_Hann_index'] - df['Mutated_Hann_index'])
df['Abs_Diff_Wiener_index'] = abs(df['Original_Wiener_index'] - df['Mutated_Wiener_index'])
df['Abs_Diff_SAS_index'] = abs(df['Original_SAS_index'] - df['Mutated_SAS_index'])
df['Abs_Diff_QED_index'] = abs(df['Original_QED_index'] - df['Mutated_QED_index'])

# Group by 'Method' and 'Mutant_Number', then calculate the average absolute difference
avg_diff = df.groupby(['Method', 'Mutant_Number']).agg({'Abs_Diff_Bertz_index': 'mean',
                                                        'Abs_Diff_Hann_index': 'mean',
                                                        'Abs_Diff_Wiener_index': 'mean',
                                                        'Abs_Diff_SAS_index': 'mean',
                                                        'Abs_Diff_QED_index': 'mean'})

## Calculate the overall average for each method
avg_diff['Average'] = avg_diff.mean(axis=1)

# Reset the index for better formatting
avg_diff.reset_index(inplace=True)

# Pivot the table to the desired format
result = avg_diff.pivot_table(index=['Method'], columns='Mutant_Number', values=['Abs_Diff_Bertz_index', 
                                                                                'Abs_Diff_Hann_index', 
                                                                                'Abs_Diff_Wiener_index', 
                                                                                'Abs_Diff_SAS_index', 
                                                                                'Abs_Diff_QED_index', 
                                                                                'Average'])

# Rename columns for better readability
result.columns = [f'{col[0]}_{col[1]}' for col in result.columns]

# Add 'Average' row to the table
result.loc['Average'] = result.mean()

# Print the result for inspection
print(result)

# Save the result to a CSV file
result.to_csv("Pairs_cleaned_complexity_indices_AveDif.csv", index=True)