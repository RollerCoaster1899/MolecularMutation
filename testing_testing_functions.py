import pandas as pd
from conversion.smiles_to_ECFP4 import smiles_to_ECFP4
from testing.validity_of_generated_molecules.validity_check import IsCorrectSMILES
from testing.conservation_relevant_features.pic50_prediction import predict_pIC50
from testing.diversity_and_complexity.molecule_complexity_indexs import bertz_index, hann_index, wiener_index, sas_index, qed_index 
from testing.diversity_and_complexity.molecule_complexity_descriptors import calculate_fraction_chiral_carbons, calculate_fraction_sp3, calculate_molecular_weight, calculate_number_of_heteroatoms, calculate_number_of_rings
from testing.diversity_and_complexity.diversity_CDF_tanimoto_ecfp4 import plot_similarity_cdf

smiles = 'CNC(C)CC1=CC=C2C(=C1)OCO2'
ecfp4 = smiles_to_ECFP4(smiles)

print("------------------------------------------------")

### Testing validity of generated molecules
# Test smiles string validity 0 = It is not valid - 1 = It is valid
validity = IsCorrectSMILES(smiles)
print("Predicted validity: ", validity)

print("------------------------------------------------")
### Testing conservation of relevant features
# Predict the pIC50 from ECFP4
predicted_pCI50 = predict_pIC50(ecfp4, "CHEMBL267_RF")
print("Predicted CHEMBL267 pIC50: ", str(predicted_pCI50))

predicted_pCI50 = predict_pIC50(ecfp4, "CHEMBL1957_RF")
print("Predicted CHEMBL1957 pIC50: ", str(predicted_pCI50))

predicted_pCI50 = predict_pIC50(ecfp4, "CHEMBL2842_RF")
print("Predicted CHEMBL2842 pIC50: ", str(predicted_pCI50))

print("------------------------------------------------")

### Testing complexity
## Testing complexity index
# Testing complexity Bertz index - Molecular topology
bertz_index = bertz_index(smiles)
print("Bertz index: ", str(bertz_index))

# Testing complexity Hann index - substructure and conectivity complexity
hann_index = hann_index(smiles)
print("Hann index: ", str(hann_index))

# Testing complexity Wiener index - molecular branching and ciclicity 
wiener_index = wiener_index(smiles)
print("Wiener index: ", str(wiener_index))

## Estimators and Scores
# Testing Quantitative Drug Estimation
sas_index = sas_index(smiles)
print("SAS index: ", str(sas_index))

# Testing Synthetic Acesbility Score
qed_index = qed_index(smiles)
print("QED index: ", str(qed_index))

## Testing complexity descriptors
# Testing Molecular Weight 
molecular_weight = calculate_molecular_weight(smiles)
print("Molecular Weight:", molecular_weight)

# Testing Number of Chiral Carbons
No_chiral_carbons = calculate_fraction_chiral_carbons(smiles)
print("Fraction of Chiral Carbons:", No_chiral_carbons)

# Testing Number of Rings
No_rings = calculate_number_of_rings(smiles)
print("Number of Rings:", No_rings)

# Testing Number of Heteroatoms
No_heteroatoms = calculate_number_of_heteroatoms(smiles)
print("Number of Heteroatoms:", No_heteroatoms)

## Combining diversity and complexity
# Testing sp3 hybridization
sp3_hybridization = calculate_fraction_sp3(smiles)
print("Fraction of sp3 Hybridization:", sp3_hybridization)

print("------------------------------------------------")

### Testing complexity
## Cumulative distributions functions
smiles_list_1 = pd.read_csv("/home/raul-acosta/GitHub/mutation_and_recombination/testing/diversity_and_complexity/comparison_datasets/fda.csv")
smiles_list_1 = smiles_list_1["smiles"]

smiles_list_2 = pd.read_csv("/home/raul-acosta/GitHub/mutation_and_recombination/testing/diversity_and_complexity/comparison_datasets/world.csv")
smiles_list_2 = smiles_list_2["smiles"]

smiles_list_3 = pd.read_csv("/home/raul-acosta/GitHub/mutation_and_recombination/testing/diversity_and_complexity/comparison_datasets/agent.csv")
smiles_list_3 = smiles_list_3["smiles"]

smiles_list_3 = pd.read_csv("/home/raul-acosta/GitHub/mutation_and_recombination/testing/diversity_and_complexity/comparison_datasets/agent.csv")
smiles_list_3 = smiles_list_3["smiles"]

smiles_list_4 = pd.read_csv("/home/raul-acosta/GitHub/mutation_and_recombination/testing/diversity_and_complexity/comparison_datasets/investigational-only.csv")
smiles_list_4 = smiles_list_4["smiles"]

labels = ["fda", "world", "agent", "investigational_only"]

plot_similarity_cdf([smiles_list_1, smiles_list_2, smiles_list_3], labels=labels)


