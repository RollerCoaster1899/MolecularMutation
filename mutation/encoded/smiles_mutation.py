# Code from Krenn, M., Florian Häse, Nigam, A., & Friederich, P. (2020). Self-referencing embedded strings (SELFIES): A 100% robust molecular string representation. Machine Learning: Science and Technology, 1(4), 045024–045024. https://doi.org/10.1088/2632-2153/aba947

import random

def mutate_smiles(input_smiles, num_mutants=1, mutation_prob=0.2, smiles_symbols=None):
    if smiles_symbols is None:
        smiles_symbols = 'FONC()=#12345'

    mutants = set()  # Using a set to store unique mutants

    while len(mutants) < num_mutants:
        mutated_smiles = input_smiles
        for _ in range(len(input_smiles)):
            if random.random() < mutation_prob:
                mol_idx = random.randint(0, len(mutated_smiles) - 1)
                symbol_idx = random.randint(0, len(smiles_symbols) - 1)
                mutated_smiles = mutated_smiles[:mol_idx] + smiles_symbols[symbol_idx] + mutated_smiles[mol_idx + 1:]
        
        # Ensure mutated_smiles is not None before adding to the set
        if mutated_smiles is not None:
            mutants.add(mutated_smiles)

    return list(mutants)  # Convert the set back to a list

