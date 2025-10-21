import numpy
from selfies import encoder, decoder
from conversion.smiles_to_ECFP4 import smiles_to_ECFP4
from conversion.selfies_to_smiles import selfies_to_smiles
from mutation.binary_vectors.bitflip_mutation import bitflip_mutator
from mutation.binary_vectors.resetting_mutation import random_resetting_mutator
from mutation.binary_vectors.swap_mutation import swap_mutator
from mutation.binary_vectors.inversion_mutation import inversion_mutator
from mutation.binary_vectors.scramble_mutation import scramble_mutator
from mutation.encoded.smiles_mutation import mutate_smiles
from mutation.encoded.selfies_mutation import mutate_selfies
from mutation.string.smilesclickchem import smiles_click_chem_mutator
from mutation.graph.gb_gm import gb_gm_mutator
from mutation.graph.gb_ga import gb_ga_mutator


# Generate morgan fingerprints
smiles = 'CN(Cc1cnc2nc(N)nc(N)c2n1)c1ccc(C(=O)N[C@@H](CCC(=O)O)C(=O)O)cc1'
selfies = encoder(smiles)
ecfp4_mf = smiles_to_ECFP4(smiles)
print("Original array: " + str(ecfp4_mf) + " Original array sum: " + str(numpy.sum(ecfp4_mf)))

probability = 1/len(ecfp4_mf)
mutants_number = 10

### Binary Vectors 0-1s
# Bitflip mutator
bf_ecfp4_mf = bitflip_mutator(ecfp4_mf, probability, mutants_number)
print("BF mutator array: " + str(bf_ecfp4_mf) + " Mutated array sum: " + str(numpy.sum(bf_ecfp4_mf)))

# Random resetting mutator
rr_ecfp4_mf = random_resetting_mutator(ecfp4_mf, probability, mutants_number)
print("RR mutator array: " + str(rr_ecfp4_mf) + " Mutated array sum: " + str(numpy.sum(rr_ecfp4_mf)))

# Swapping mutator
sw_ecfp4_mf = swap_mutator(ecfp4_mf, probability, mutants_number)
print("SW mutator array: " + str(sw_ecfp4_mf) + " Mutated array sum: " + str(numpy.sum(sw_ecfp4_mf)))

# Scrambling mutator
sc_ecfp4_mf = scramble_mutator(ecfp4_mf, probability, mutants_number)
print("SC mutator array: " + str(sc_ecfp4_mf) + " Mutated array sum: " + str(numpy.sum(sc_ecfp4_mf)))

# Inversion mutator
im_ecfp4_mf = inversion_mutator(ecfp4_mf, probability, mutants_number)
print("IN mutator array: " + str(im_ecfp4_mf) + " Mutated array sum: " + str(numpy.sum(im_ecfp4_mf)))

### Encoded
# SMILES Mutation 
mutated_smiles = mutate_smiles(smiles, mutants_number, probability, smiles_symbols=None)
print("SM Mutator SMILES:", mutated_smiles)

# SELFIES Mutation
mutated_selfies = mutate_selfies(selfies, mutants_number, probability, selfies_symbols=None)
print("SM Mutator SELFIES:", mutated_selfies)
mutated_selfies = selfies_to_smiles(mutated_selfies)
print("SM Mutator SELFIES:", mutated_selfies)

## String
# SMILES Mutation
scc_smiles = smiles_click_chem_mutator(smiles, mutants_number)
print("SCC mutator SMILES: ", scc_smiles)

## Graph 
# GB_GM 
from rdkit import Chem

#input_data = [{"SMILES": smiles, "Num_Mutants": 3, "Max_Atoms": 15, "Average_Size": 10, "Size_Stdev": 2}]

original_molecule = Chem.MolFromSmiles(smiles)
num_atoms = original_molecule.GetNumAtoms()
# Generate mutants using the specified properties
input_data = [{"SMILES": smiles, "Num_Mutants": mutants_number, "Max_Atoms": (num_atoms*1.3), "Average_Size": num_atoms, "Size_Stdev": (num_atoms*0.1913)}]
mutated_smiles = gb_gm_mutator(input_data)
print("GB GM mutator SMILES: ", mutated_smiles)

# GB_GA
input_data = [{"SMILES": smiles, "Num_Mutants": mutants_number}]
mutated_smiles = gb_ga_mutator(input_data)
print("GB GA mutator SMILES: ", mutated_smiles)


