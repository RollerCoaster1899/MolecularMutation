# Code from Krenn, M., Florian Häse, Nigam, A., & Friederich, P. (2020). Self-referencing embedded strings (SELFIES): A 100% robust molecular string representation. Machine Learning: Science and Technology, 1(4), 045024–045024. https://doi.org/10.1088/2632-2153/aba947

from rdkit.Chem import MolFromSmiles
from selfies import encoder, decoder
import random

def tokenize_selfies(selfies):
    location = selfies.find(']')
    all_tokens = []
    while location >= 0:
        all_tokens.append(selfies[0:location + 1])
        selfies = selfies[location + 1:]
        location = selfies.find(']')
    return all_tokens

def detokenize_selfies(selfies_list):
    return ''.join(selfies_list)

# Define the function mutate_selfies
def mutate_selfies(input_selfies, num_mutants=1, mutation_prob=0.2, selfies_symbols=None):
    if selfies_symbols is None:
        selfies_symbols = ['[epsilon]', '[Ring1]', '[Ring2]', '[Branch1_1]', '[Branch1_2]', '[Branch1_3]',
                            '[F]', '[O]', '[=O]', '[N]', '[=N]', '[#N]', '[C]', '[=C]', '[#C]']

    mutants = set()  # Using a set to store unique mutants

    if isinstance(input_selfies, list):
        input_selfies = detokenize_selfies(input_selfies)

    # Loop until num_mutants unique mutants are generated
    while len(mutants) < num_mutants:
        input_selfies_tok = tokenize_selfies(input_selfies)

        for _ in range(num_mutants):
            for idx in range(len(input_selfies_tok)):
                if random.random() < mutation_prob:
                    symbol_idx = random.randint(0, len(selfies_symbols) - 1)
                    new_selfies_str = detokenize_selfies(input_selfies_tok[0:idx])
                    new_selfies_str += selfies_symbols[symbol_idx]
                    new_selfies_str += detokenize_selfies(input_selfies_tok[idx + 1:])
                    input_selfies_tok = tokenize_selfies(new_selfies_str)

            mutated_selfies = detokenize_selfies(input_selfies_tok)
            
            # Ensure mutated_selfies is not None before adding to the set
            if mutated_selfies is not None:
                mutants.add(mutated_selfies)
            
            # Break the loop if num_mutants unique mutants are generated
            if len(mutants) == num_mutants:
                break

    return list(mutants)  # Convert the set back to a list