import os
import pandas as pd
import pandas as pd
from conversion.smiles_to_ECFP4 import smiles_to_ECFP4
from testing.conservation_relevant_features.pic50_prediction import predict_pIC50
import pandas as pd
from conversion.smiles_to_ECFP4 import smiles_to_ECFP4
from testing.validity_of_generated_molecules.validity_check import IsCorrectSMILES
from testing.conservation_relevant_features.pic50_prediction import predict_pIC50
from testing.diversity_and_complexity.molecule_complexity_indexs import bertz_index, hann_index, wiener_index, sas_index, qed_index 
from testing.diversity_and_complexity.molecule_complexity_descriptors import calculate_fraction_chiral_carbons, calculate_fraction_sp3, calculate_molecular_weight, calculate_number_of_heteroatoms, calculate_number_of_rings
from testing.diversity_and_complexity.diversity_CDF_tanimoto_ecfp4 import plot_similarity_cdf


folder_path = "/home/raul-acosta/GitHub/mutation_and_recombination/results/computational_resources/Pairs.csv"
df = pd.read_csv(folder_path)

# Function to calculate complexity descriptors for each SMILES
def calculate_complexity_descriptors(smiles):
    try:
        molecular_weight = calculate_molecular_weight(smiles)  # Assuming you have defined this function
        No_chiral_carbons = calculate_fraction_chiral_carbons(smiles)  # Assuming you have defined this function
        No_rings = calculate_number_of_rings(smiles)  # Assuming you have defined this function
        No_heteroatoms = calculate_number_of_heteroatoms(smiles)  # Assuming you have defined this function
        sp3_hybridization = calculate_fraction_sp3(smiles)  # Assuming you have defined this function
        return molecular_weight, No_chiral_carbons, No_rings, No_heteroatoms, sp3_hybridization
    except ValueError as e:
        print(f"Error calculating complexity descriptors for SMILES '{smiles}': {e}")
        return None, None, None, None, None

"""

# Apply the function to calculate complexity descriptors for each Original_SMILES
df[['Original_Molecular_Weight', 'Original_Fraction_of_Chiral_Carbons', 
    'Original_Number_of_Rings', 'Original_Number_of_Heteroatoms', 
    'Original_Fraction_of_sp3_Hybridization']] = pd.DataFrame(df['Original_SMILES'].apply(calculate_complexity_descriptors).tolist(), index=df.index)


# Split the "Mutated_SMILES" column into separate rows
df = df.assign(Mutated_SMILES=df['Mutated_SMILES'].str.strip("[]").str.split(", ")).explode('Mutated_SMILES')

# Remove single quotes from each mutated SMILES
df['Mutated_SMILES'] = df['Mutated_SMILES'].str.strip("'")

# Move the "Mutated_SMILES" column to the last position
mutated_smiles_column = df.pop('Mutated_SMILES')
df['Mutated_SMILES'] = mutated_smiles_column

# Apply the function to calculate complexity descriptors for each Mutated_SMILES
df[['Mutated_Molecular_Weight', 'Mutated_Fraction_of_Chiral_Carbons', 
    'Mutated_Number_of_Rings', 'Mutated_Number_of_Heteroatoms', 
    'Mutated_Fraction_of_sp3_Hybridization']] = pd.DataFrame(df['Mutated_SMILES'].apply(calculate_complexity_descriptors).tolist(), index=df.index)
# Specify the new columns for cleaning
new_columns_to_check = ['Mutated_SMILES', 'Mutated_Molecular_Weight', 'Mutated_Fraction_of_Chiral_Carbons', 
                        'Mutated_Number_of_Rings', 'Mutated_Number_of_Heteroatoms', 
                        'Mutated_Fraction_of_sp3_Hybridization']

# Drop rows with NaN or None values in the new specified columns
df = df.dropna(subset=new_columns_to_check)

# Reset the index after dropping rows
df.reset_index(drop=True, inplace=True)

# Save the cleaned DataFrame to a new CSV file
df.to_csv("Pairs_cleaned_complexity_descriptors.csv", index=False)

"""


folder_path = "Pairs_cleaned_complexity_descriptors.csv"
df = pd.read_csv(folder_path)

# Calculate the absolute differences for each complexity descriptor
df['Abs_Diff_Molecular_Weight'] = abs(df['Original_Molecular_Weight'] - df['Mutated_Molecular_Weight'])
df['Abs_Diff_Fraction_of_Chiral_Carbons'] = abs(df['Original_Fraction_of_Chiral_Carbons'] - df['Mutated_Fraction_of_Chiral_Carbons'])
df['Abs_Diff_Number_of_Rings'] = abs(df['Original_Number_of_Rings'] - df['Mutated_Number_of_Rings'])
df['Abs_Diff_Number_of_Heteroatoms'] = abs(df['Original_Number_of_Heteroatoms'] - df['Mutated_Number_of_Heteroatoms'])
df['Abs_Diff_Fraction_of_sp3_Hybridization'] = abs(df['Original_Fraction_of_sp3_Hybridization'] - df['Mutated_Fraction_of_sp3_Hybridization'])

# Save the DataFrame with absolute differences to a new CSV file
df.to_csv("Pairs_cleaned_complexity_descriptors_difference.csv", index=False)

# Calculate the absolute differences for each complexity descriptor
df['Abs_Diff_Molecular_Weight'] = abs(df['Original_Molecular_Weight'] - df['Mutated_Molecular_Weight'])
df['Abs_Diff_Fraction_of_Chiral_Carbons'] = abs(df['Original_Fraction_of_Chiral_Carbons'] - df['Mutated_Fraction_of_Chiral_Carbons'])
df['Abs_Diff_Number_of_Rings'] = abs(df['Original_Number_of_Rings'] - df['Mutated_Number_of_Rings'])
df['Abs_Diff_Number_of_Heteroatoms'] = abs(df['Original_Number_of_Heteroatoms'] - df['Mutated_Number_of_Heteroatoms'])
df['Abs_Diff_Fraction_of_sp3_Hybridization'] = abs(df['Original_Fraction_of_sp3_Hybridization'] - df['Mutated_Fraction_of_sp3_Hybridization'])

# Group by 'Method' and 'Mutant_Number', then calculate the average absolute difference
avg_diff = df.groupby(['Method', 'Mutant_Number']).agg({'Abs_Diff_Molecular_Weight': 'mean',
                                                        'Abs_Diff_Fraction_of_Chiral_Carbons': 'mean',
                                                        'Abs_Diff_Number_of_Rings': 'mean',
                                                        'Abs_Diff_Number_of_Heteroatoms': 'mean',
                                                        'Abs_Diff_Fraction_of_sp3_Hybridization': 'mean'})

# Calculate the overall average for each method
avg_diff['Average'] = avg_diff.mean(axis=1)

# Reset the index for better formatting
avg_diff.reset_index(inplace=True)

# Pivot the table to the desired format
result = avg_diff.pivot_table(index='Method', columns='Mutant_Number', values=['Abs_Diff_Molecular_Weight', 
                                                                                'Abs_Diff_Fraction_of_Chiral_Carbons', 
                                                                                'Abs_Diff_Number_of_Rings', 
                                                                                'Abs_Diff_Number_of_Heteroatoms', 
                                                                                'Abs_Diff_Fraction_of_sp3_Hybridization', 
                                                                                'Average'])

# Rename columns for better readability
result.columns = [f'{col[0]}_{col[1]}' for col in result.columns]

# Add 'Average' row to the table
result.loc['Average'] = result.mean()

# Print the result for inspection
print(result)

# Save the result to a CSV file
result.to_csv("Pairs_cleaned_complexity_descriptors_AveDif.csv")
