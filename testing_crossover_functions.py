import numpy
from selfies import encoder, decoder
from conversion.smiles_to_ECFP4 import smiles_to_ECFP4
from crossover.binary_vectors.one_point_crossover import one_point_crossover
from crossover.binary_vectors.multi_point_crossover import multi_point_crossover
from crossover.binary_vectors.uniform_crossover import uniform_crossover
from crossover.binary_vectors.whole_arithmetic_crossover import whole_arithmetic_recombination
from crossover.binary_vectors.partially_mapped_crossover import partially_mapped_crossover
from crossover.binary_vectors.order_based_crossover import order_based_crossover
from crossover.binary_vectors.shuffle_crossover import shuffle_crossover
from crossover.binary_vectors.ring_crossover import ring_crossover
from crossover.string.smilesmerge import smiles_merge_crossover
from crossover.string.gb_ga import gb_ga_crossover

parent1_smiles = "CCc1ccc(N=[N+]=[N-])cc1"
parent2_smiles = "O=C(C)Oc1ccccc1C(=O)O"
parents_smiles = [parent1_smiles, parent2_smiles]

ecfp4_mf_p1 = smiles_to_ECFP4(parent1_smiles)
print("Original array: " + str(ecfp4_mf_p1) + " Original array sum: " + str(numpy.sum(ecfp4_mf_p1)))

ecfp4_mf_p2 = smiles_to_ECFP4(parent2_smiles)
print("Original array: " + str(ecfp4_mf_p2) + " Original array sum: " + str(numpy.sum(ecfp4_mf_p2)))

parent1 = ecfp4_mf_p1
parent2 = ecfp4_mf_p2

crossover_probability = 0.8
crossover_point_probability = 1/len(ecfp4_mf_p1)
bias = None #0.6 Optional bias towards parent1
alpha = 0.5  # Weighting factor, use 0.5 for identical children
num_offspring = 1

# One point crossover 
op_ecpf_mf = one_point_crossover(parent1, parent2, crossover_probability, crossover_point_probability, num_offspring)
print("op point crossover array: " + str(op_ecpf_mf))

# Multi point crossover
mp_ecpf_mf = multi_point_crossover(parent1, parent2, crossover_probability, num_offspring)
print("Multi point crossover array: " + str(mp_ecpf_mf))

# Uniform crossover
uc_ecpf_mf = uniform_crossover(parent1, parent2, crossover_probability, bias, num_offspring)
print("Uniform crossover array: " + str(uc_ecpf_mf))

# Whole aritmetic crossover # don´t work for mf as obtains the average 
wa_ecpf_mf = whole_arithmetic_recombination(parent1, parent2, alpha, num_offspring)
print("Whole aritmetic crossover array: " + str(wa_ecpf_mf))

# Partially mapped crossover
pm_ecpf_mf = partially_mapped_crossover(parent1, parent2, crossover_probability, num_offspring)
print("Partially mapped crossover array: " + str(pm_ecpf_mf))

# Shuffle Crossover
sc_ecpf_mf = shuffle_crossover(parent1, parent2, crossover_probability, num_offspring)
print("Partially mapped crossover array: " + str(sc_ecpf_mf))

# Ring crossover
rc_ecpf_mf = ring_crossover(parent1, parent2, crossover_probability, num_offspring)
print("Ring crossover array: " + str(rc_ecpf_mf))

## String
# SMILESMerge crossover
sm_smiles = smiles_merge_crossover(parents_smiles, num_offspring)
print("SmilesMerge crossover: ", sm_smiles)

## GB-GA crossover
input_data = [{"SMILES_1": parent1_smiles, "SMILES_2": parent2_smiles, "Num_children": 1}]
children = gb_ga_crossover(input_data)
print(children)