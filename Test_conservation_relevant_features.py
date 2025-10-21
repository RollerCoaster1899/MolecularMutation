import os
import pandas as pd
import pandas as pd
from conversion.smiles_to_ECFP4 import smiles_to_ECFP4
from testing.conservation_relevant_features.pic50_prediction import predict_pIC50


folder_path = "/home/raul-acosta/GitHub/mutation_and_recombination/results/computational_resources/Pairs.csv"
df = pd.read_csv(folder_path)
#df = df.tail(5)

# Function to calculate pIC50 for each SMILES using a specified model
def calculate_pIC50(smiles, index, model_name):
    # Check if the SMILES is enclosed in square brackets and single quotes
    if smiles.startswith("['") and smiles.endswith("']"):
        mutated_smiles = smiles[2:-2].split("', '")  # Extract individual mutated SMILES
        pIC50_values = []
        for mutated_smile in mutated_smiles:
            ecfp4 = smiles_to_ECFP4(mutated_smile)
            if ecfp4 is not None:
                print(f"Predicting pIC50 for mutated SMILES '{mutated_smile}' at index {index} using model {model_name}")
                pIC50_value = predict_pIC50(ecfp4, model_name)  # Assuming you have defined this function
                if pIC50_value is not None:
                    pIC50_values.append(pIC50_value)
                else:
                    print(f"Error predicting pIC50 for mutated SMILES '{mutated_smile}' at index {index}")
                    pIC50_values.append(None)
            else:
                print(f"Error: Failed to parse mutated SMILES '{mutated_smile}' at index {index}")
                pIC50_values.append(None)
        return pIC50_values
    else:
        ecfp4 = smiles_to_ECFP4(smiles)
        if ecfp4 is not None:
            print(f"Predicting pIC50 for SMILES at index {index} using model {model_name}")
            pIC50_value = predict_pIC50(ecfp4, model_name)  # Assuming you have defined this function
            if pIC50_value is not None:
                return pIC50_value
            else:
                print(f"Error predicting pIC50 for SMILES '{smiles}' at index {index}")
                return None
        else:
            print(f"Error: Failed to parse SMILES '{smiles}' at index {index}")
            return None
"""

# Apply the function to calculate pIC50 for each Original_SMILES using different models
df['Original_pIC50_CHEMBL267'] = df.apply(lambda row: calculate_pIC50(row['Original_SMILES'], row.name, "CHEMBL267_RF"), axis=1)
df['Original_pIC50_CHEMBL1957'] = df.apply(lambda row: calculate_pIC50(row['Original_SMILES'], row.name, "CHEMBL1957_RF"), axis=1)
df['Original_pIC50_CHEMBL2842'] = df.apply(lambda row: calculate_pIC50(row['Original_SMILES'], row.name, "CHEMBL2842_RF"), axis=1)


# Split the "Mutated_SMILES" column into separate rows
df = df.assign(Mutated_SMILES=df['Mutated_SMILES'].str.strip("[]").str.split(", ")).explode('Mutated_SMILES')

# Remove single quotes from each mutated SMILES
df['Mutated_SMILES'] = df['Mutated_SMILES'].str.strip("'")

# Move the "Mutated_SMILES" column to the last position
mutated_smiles_column = df.pop('Mutated_SMILES')
df['Mutated_SMILES'] = mutated_smiles_column


# Predict pIC50 for the mutated SMILES and append the predictions to the DataFrame
df['Mutated_pIC50_CHEMBL267'] = df.apply(lambda row: calculate_pIC50(row['Mutated_SMILES'], row.name, "CHEMBL267_RF"), axis=1)
df['Mutated_pIC50_CHEMBL1957'] = df.apply(lambda row: calculate_pIC50(row['Mutated_SMILES'], row.name, "CHEMBL1957_RF"), axis=1)
df['Mutated_pIC50_CHEMBL2842'] = df.apply(lambda row: calculate_pIC50(row['Mutated_SMILES'], row.name, "CHEMBL2842_RF"), axis=1)


"""


folder_path = "/home/raul-acosta/GitHub/mutation_and_recombination/Pairs_updated3.csv"
df = pd.read_csv(folder_path)

# Drop rows with NaN or None values in the specified columns
columns_to_check = ['Mutated_SMILES', 'Mutated_pIC50_CHEMBL267', 'Mutated_pIC50_CHEMBL1957', 'Mutated_pIC50_CHEMBL2842']
df = df.dropna(subset=columns_to_check)

# Reset the index after dropping rows
df.reset_index(drop=True, inplace=True)

df.to_csv("Pairs_cleaned_pIC50.csv")

# Calculate the absolute differences for each model
df['Abs_Diff_CHEMBL267'] = abs(df['Original_pIC50_CHEMBL267'] - df['Mutated_pIC50_CHEMBL267'])
df['Abs_Diff_CHEMBL1957'] = abs(df['Original_pIC50_CHEMBL1957'] - df['Mutated_pIC50_CHEMBL1957'])
df['Abs_Diff_CHEMBL2842'] = abs(df['Original_pIC50_CHEMBL2842'] - df['Mutated_pIC50_CHEMBL2842'])

# Save the DataFrame to a new CSV file
df.to_csv("Pairs_cleaned_pIC50_difference.csv", index=False)

# Calculate the absolute differences for each method
df['Abs_Diff_CHEMBL267'] = abs(df['Original_pIC50_CHEMBL267'] - df['Mutated_pIC50_CHEMBL267'])
df['Abs_Diff_CHEMBL1957'] = abs(df['Original_pIC50_CHEMBL1957'] - df['Mutated_pIC50_CHEMBL1957'])
df['Abs_Diff_CHEMBL2842'] = abs(df['Original_pIC50_CHEMBL2842'] - df['Mutated_pIC50_CHEMBL2842'])

# Group by Method and number of mutants, then calculate the average absolute difference
avg_diff = df.groupby(['Method', 'Mutant_Number']).agg({'Abs_Diff_CHEMBL267': 'mean', 
                                                        'Abs_Diff_CHEMBL1957': 'mean', 
                                                        'Abs_Diff_CHEMBL2842': 'mean'})

# Calculate the overall average for each method
avg_diff['Average'] = avg_diff.mean(axis=1)

# Reset the index for better formatting
avg_diff.reset_index(inplace=True)

# Pivot the table to the desired format
result = avg_diff.pivot_table(index='Method', columns='Mutant_Number', values=['Abs_Diff_CHEMBL267', 'Abs_Diff_CHEMBL1957', 'Abs_Diff_CHEMBL2842', 'Average'])

# Rename columns for better readability
result.columns = [f'{col[0]}_{col[1]}' for col in result.columns]

# Add 'Average' row to the table
result.loc['Average'] = result.mean()

print(result)

result.to_csv("Pairs_cleaned_pIC50_AveDif.csv")