import pandas as pd
from testing.diversity_and_complexity.diversity_CDF_tanimoto_ecfp4 import plot_similarity_cdf
import matplotlib.pyplot as plt
from tqdm import tqdm

# Read the CSV file
pairs_file = "/home/raul-acosta/GitHub/mutation_and_recombination/results/conservation_relevant_features/Pairs_cleaned_pIC50.csv"
df = pd.read_csv(pairs_file)

# Group by 'Method' and collect the original and mutant SMILES for each method
method_smiles_dict = {}
for method, group in tqdm(df.groupby('Method'), desc="Processing Methods", leave=False):
    smiles_lists = [group[group['Mutant_Number'] == f'Mutant_{i}']['Original_SMILES'].tolist() for i in [1, 3, 5]]
    smiles_lists.append(group['Mutated_SMILES'].tolist())  # Add the Mutated_SMILES for the mutants
    method_smiles_dict[method] = smiles_lists

# Define labels
labels = ["Original", "Mutant 1", "Mutant 3", "Mutant 5", "Mutated"]

# Plot and save each plot for each method
for method, smiles_lists in tqdm(method_smiles_dict.items(), desc="Processing Methods", leave=False):
    plot_similarity_cdf(smiles_lists, labels=labels, title=f'Similarity Distribution for Method: {method}')
    plt.savefig(f"Similarity_CDF_{method}.png")
    plt.close()

# Plot all methods together in one image
fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(15, 10))
for (method, smiles_lists), ax in tqdm(zip(method_smiles_dict.items(), axes.flatten()), desc="Processing Methods", leave=False):
    plot_similarity_cdf(smiles_lists, labels=labels, title=f'Similarity Distribution for Method: {method}', ax=ax)

plt.tight_layout()
plt.savefig("Similarity_CDF_All_Methods.png")
plt.close()
