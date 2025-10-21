import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def calculate_tanimoto_similarity(smiles_list_1, smiles_list_2):
    mols_1 = [Chem.MolFromSmiles(smiles) for smiles in smiles_list_1]
    mols_2 = [Chem.MolFromSmiles(smiles) for smiles in smiles_list_2]
    fingerprints_1 = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048) for mol in mols_1]
    fingerprints_2 = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048) for mol in mols_2]

    # Calculate pairwise Tanimoto coefficients
    tanimoto_similarities = []
    for fp_1 in fingerprints_1:
        for fp_2 in fingerprints_2:
            tanimoto_similarity = AllChem.DataStructs.TanimotoSimilarity(fp_1, fp_2)
            tanimoto_similarities.append(tanimoto_similarity)

    return np.array(tanimoto_similarities)

def plot_similarity_violin(smiles_lists, labels=None, title=None):
    # Calculate Tanimoto similarities for each set of SMILES
    similarities = []
    for smiles_list in smiles_lists:
        original_smiles = smiles_list[0]
        mutated_smiles = smiles_list[-1]
        similarity = calculate_tanimoto_similarity(original_smiles, mutated_smiles)
        similarities.extend(similarity)

    # Create a DataFrame for plotting
    df = pd.DataFrame(similarities, columns=["Tanimoto Similarity"])
    df["Method"] = labels[0]

    # Plotting
    plt.figure(figsize=(10, 6))
    sns.violinplot(x="Method", y="Tanimoto Similarity", data=df, palette="viridis")

    # Customize plot
    plt.xlabel("Method")
    plt.ylabel("Tanimoto Similarity")
    if title:
        plt.title(title)
    plt.savefig("Similarity_Violin_Plot.png", dpi=300)
    plt.show()

# Example usage:
# Read the data
pairs_file = "/home/raul-acosta/GitHub/mutation_and_recombination/results/conservation_relevant_features/Pairs_cleaned_pIC50.csv"
df = pd.read_csv(pairs_file)

# Group by 'Method' and collect the original and mutant SMILES for each method
method_smiles_dict = {}
for method, group in df.groupby('Method'):
    original_smiles = group[group['Mutant_Number'] == 'Original']['Original_SMILES'].tolist()
    mutated_smiles = group[group['Mutant_Number'] == 'Mutated']['Mutated_SMILES'].tolist()
    method_smiles_dict[method] = [original_smiles, mutated_smiles]

# Define labels
labels = list(method_smiles_dict.keys())

# Plot and save the violin plot
plot_similarity_violin(method_smiles_dict.values(), labels=labels)
