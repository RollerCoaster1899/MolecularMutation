import pandas as pd
from time import time
from rdkit import Chem
from tqdm import tqdm
from testing.validity_of_generated_molecules.validity_check import IsCorrectSMILES
from mutation.encoded.smiles_mutation import mutate_smiles
from conversion.selfies_to_smiles import selfies_to_smiles
from conversion.smiles_to_selfies import smiles_to_selfies
from mutation.encoded.selfies_mutation import mutate_selfies
import numpy
import random
from mutation.encoded.smiles_mutation import mutate_smiles
from mutation.encoded.selfies_mutation import mutate_selfies
from mutation.string.smilesclickchem import smiles_click_chem_mutator
from mutation.graph.gb_gm import gb_gm_mutator
from itertools import chain
import multiprocessing
from multiprocessing import TimeoutError, Manager, Process
from mutation.graph.gb_ga import gb_ga_mutator
import os



# Open the file and read the SMILES strings into a list



results_path = "/home/raul-acosta/GitHub/mutation_and_recombination/results/computational_resources"

# Open the file and read the SMILES strings into a list
file_path = "/home/raul-acosta/GitHub/mutation_and_recombination/results/Random_initial_sample.txt"
with open(file_path, 'r') as file:
    smiles_list = file.readlines()

# Strip newline characters from each SMILES string
smiles_list = [smiles.strip() for smiles in smiles_list]

# Now you have the list of SMILES strings
print(smiles_list)



"""

df = pd.read_csv("/home/raul-acosta/GitHub/mutation_and_recombination/testing/diversity_and_complexity/comparison_datasets/fda.csv")
smiles_list = random.sample(list(df["smiles"]), 500)

# Now you have the list of SMILES strings
print(smiles_list)



# Open a text file in write mode
filename = f'Random_initial_sample.txt'
# Join the path and filename
file_path = os.path.join(results_path, filename)
# Open the file for writing
with open(file_path, 'w') as f:
    # Write each mutant to the file
    for sample in smiles_list:
        f.write(sample + '\n')   

        
"""

# ---------------------------------------------
## SMILES 

# Initialize lists to store results
time_taken = []
percent_validity = []

# Initialize headers for CSV files
headers = [1, 3, 5]

# Iterate over each mutant number
for mutants_number in tqdm(headers, desc="Generating mutants"):
    # Start overall time measurement
    overall_start_time = time()
    all_mutants = []
    pairs = []

    # Generate mutants for all SMILES in the list
    for smiles in smiles_list:
        try:
            probability = 1 / len(smiles)

            # Generate mutants using the specified properties
            mutated_smiles = mutate_smiles(smiles, mutants_number, probability, smiles_symbols=None)
            all_mutants.extend(mutated_smiles)

            # Check if mutants are generated
            if mutated_smiles:
                # Append the pair of original SMILES and its mutants to the list
                pairs.append((smiles, mutated_smiles))

        except Exception as e:
            print(f"Error processing SMILES: {smiles}")
            print(f"Error details: {e}")
            continue

    # End overall time measurement
    overall_end_time = time()
    overall_time_taken = overall_end_time - overall_start_time
    time_taken.append(overall_time_taken)

    # Calculate percentual validity for all mutants generated with the current mutant number
    valid_count = sum(IsCorrectSMILES(mutant) for mutant in all_mutants)
    percent_validity.append(valid_count / (mutants_number * len(smiles_list)) * 100)

    # Open a text file in write mode for writing mutants
    filename = f'smiles_mutants_{mutants_number}.txt'
    file_path = os.path.join(results_path, filename)
    with open(file_path, 'w') as f:
        for mutant in all_mutants:
            f.write(mutant + '\n')

    # Open a text file in write mode for writing pairs
    filename = f"smiles_pairs_{mutants_number}.txt"
    pairs_file_path = os.path.join(results_path, filename)
    with open(pairs_file_path, "w") as file:
        for pair in pairs:
            original_smiles, mutants = pair
            file.write("Original SMILES: " + original_smiles + "\n")
            file.write("Mutants:\n")
            for mutant in mutants:
                if mutant is not None:  # Check if mutant is not None
                    file.write(mutant + "\n")
            file.write("\n")  # Add a blank line between pairs


# Create DataFrames and save as CSVs (transposed) with specified headers
time_df = pd.DataFrame({"Time Taken (seconds)": time_taken}, index=headers).transpose()
validity_df = pd.DataFrame({"Percentual Validity": percent_validity}, index=headers).transpose()

# Define the filenames for CSVs
time_filename = "smiles_time_taken.csv"
validity_filename = "smiles_percentual_validity.csv"

# Define the file paths for CSVs
time_file_path = os.path.join(results_path, time_filename)
validity_file_path = os.path.join(results_path, validity_filename)

# Save DataFrames to CSVs with specified paths and filenames
time_df.to_csv(time_file_path, index=False)
validity_df.to_csv(validity_file_path, index=False)

# Print overall time taken
print("Overall Time Taken (seconds):", overall_time_taken)

# ---------------------------------------------
## Initialize lists to store results
time_taken = []
percent_validity = []

# Initialize headers for CSV files
headers = [1, 3, 5]

# Define a function to process each molecule
def process_molecule(smiles, mutants_number, time_limit, output_list, pairs):
    try:
        probability = 1 / len(smiles)

        # Encode using SELFIES
        selfies = smiles_to_selfies(smiles)
        # Mutate using SELFIES
        mutated_selfies = mutate_selfies(selfies, mutants_number, probability, selfies_symbols=None)
        # Convert SELFIES back to SMILES
        smiles_multiple = selfies_to_smiles(mutated_selfies)

        # Filter out None values before adding to the list
        output_list.extend(filter(None, smiles_multiple))

        # If processing completed successfully, append the pair of original SMILES and its mutants to the list
        pairs.append((smiles, smiles_multiple))

    except Exception as e:
        print(f"Error processing SMILES: {smiles}")
        print(f"Error details: {e}")

# Iterate over each mutant number
for mutants_number in tqdm(headers, desc="Generating mutants"):
    # Start overall time measurement
    overall_start_time = time()
    all_mutants = []
    pairs = []

    # Create a manager to share lists among processes
    with Manager() as manager:
        output_list = manager.list()
        pair_list = manager.list()

        # Create a list to store each process
        processes = []

        i = 0
        # Process each molecule in a separate process
        for smiles in smiles_list:
            i += 1
            print("Method: SELFIES Index: ", i)
            
            process = Process(target=process_molecule, args=(smiles, mutants_number, 15, output_list, pair_list))
            process.start()
            processes.append(process)

        # Wait for each process to finish or timeout
        for process in processes:
            process.join(timeout=15)
            if process.is_alive():
                process.terminate()
                print("Timeout: Process terminated.")

        # Collect results from the shared lists
        all_mutants.extend(output_list)
        pairs.extend(pair_list)

    # End overall time measurement
    overall_end_time = time()
    overall_time_taken = overall_end_time - overall_start_time
    time_taken.append(overall_time_taken)

    # Calculate percentual validity for all mutants generated with the current mutant number
    valid_count = sum(IsCorrectSMILES(mutant) for mutant in all_mutants)
    percent_validity.append(valid_count / (mutants_number * len(smiles_list)) * 100)

    # Open a text file in write mode
    filename = f'selfies_mutants_{mutants_number}.txt'
    file_path = os.path.join(results_path, filename)
    with open(file_path, 'w') as f:
        # Write each mutant to the file
        for mutant in all_mutants:
            f.write(mutant + '\n')

# Open a text file in write mode for pairs
    filename = f'selfies_pairs_{mutants_number}.txt'
    file_path = os.path.join(results_path, filename)
    with open(file_path, "w") as file:  # Use file_path instead of pairs_file_path
        for pair in pairs:
            original_smiles, mutants = pair
            file.write("Original SMILES: " + original_smiles + "\n")
            file.write("Mutants:\n")
            for mutant in mutants:
                if mutant is not None:  # Check if mutant is not None
                    file.write(mutant + "\n")
            file.write("\n")  # Add a blank line between pairs


# Create DataFrames and save as CSVs (transposed) with specified headers
time_df = pd.DataFrame({"Time Taken (seconds)": time_taken}, index=headers).transpose()
validity_df = pd.DataFrame({"Percentual Validity": percent_validity}, index=headers).transpose()

# Define the filenames
time_filename = "selfies_time_taken.csv"
validity_filename = "selfies_percentual_validity.csv"

# Join the path and filenames
time_file_path = os.path.join(results_path, time_filename)
validity_file_path = os.path.join(results_path, validity_filename)

# Save time_df to CSV with specified path and filename
time_df.to_csv(time_file_path, index=False)

# Save validity_df to CSV with specified path and filename
validity_df.to_csv(validity_file_path, index=False)

# ---------------------------------------------

## SCC

# Define headers for CSV files
headers = [1, 3, 5]

def process_molecule(smiles, mutants_number, all_mutants, lock, pairs):
    try:
        # Generate mutants using the specified properties
        scc_smiles = smiles_click_chem_mutator(smiles, mutants_number)
        
        with lock:
            all_mutants.extend(scc_smiles)
            # Check if mutants are generated
            if scc_smiles:
                # Append the pair of original SMILES and its mutants to the list
                pairs.append((smiles, scc_smiles))
    except Exception as e:
        print(f"Error processing SMILES: {smiles}")
        print(f"Error details: {e}")

# Initialize lists to store results
time_taken = []
percent_validity = []
all_pairs = []

# Iterate over each mutant number
for mutants_number in tqdm(headers, desc="Generating mutants"):
    # Start overall time measurement
    overall_start_time = time()

    # Initialize Manager and Lock for shared resources
    all_mutants = Manager().list()
    lock = Manager().Lock()
    pairs = Manager().list()

    # Generate mutants for all SMILES in the list using multiprocessing
    processes = []
    for smiles in smiles_list:
        process = Process(target=process_molecule, args=(smiles, mutants_number, all_mutants, lock, pairs))
        processes.append(process)
        process.start()

        # Wait for the process to finish or timeout (15 seconds)
        process.join(timeout=15)
        if process.is_alive():
            process.terminate()

    # End overall time measurement
    overall_end_time = time()
    overall_time_taken = overall_end_time - overall_start_time
    time_taken.append(overall_time_taken)

    # Calculate percentual validity for all mutants generated with the current mutant number
    valid_count = sum(IsCorrectSMILES(mutant) for mutant in all_mutants)
    percent_validity.append(valid_count / (mutants_number * len(smiles_list)) * 100)

    # Save mutants to file
    mutants_filename = f'SCC_mutants_{mutants_number}.txt'
    mutants_file_path = os.path.join(results_path, mutants_filename)
    with open(mutants_file_path, 'w') as f:
        for mutant in all_mutants:
            f.write(mutant + '\n')

    # Append pairs to all_pairs list
    all_pairs.append(pairs)

# Save pairs to files
    pairs_filename = f"SCC_pairs_{mutants_number}.txt"
    pairs_file_path = os.path.join(results_path, pairs_filename)
    with open(pairs_file_path, "w") as file:
        for pair in pairs:
            original_smiles, mutants = pair
            file.write("Original SMILES: " + original_smiles + "\n")
            file.write("Mutants:\n")
            for mutant in mutants:
                file.write(mutant + "\n")
            file.write("\n")  # Add a blank line between pairs


# Create DataFrames and save as CSVs (transposed) with specified headers
time_df = pd.DataFrame({"Time Taken (seconds)": time_taken}, index=headers).transpose()
validity_df = pd.DataFrame({"Percentual Validity": percent_validity}, index=headers).transpose()

# Define the filenames
time_filename = "SCC_time_taken.csv"
validity_filename = "SCC_percentual_validity.csv"

# Save time_df to CSV with specified path and filename
time_file_path = os.path.join(results_path, time_filename)
time_df.to_csv(time_file_path, index=False)

# Save validity_df to CSV with specified path and filename
validity_file_path = os.path.join(results_path, validity_filename)
validity_df.to_csv(validity_file_path, index=False)

# Print overall time taken
print("Overall Time Taken (seconds):", overall_time_taken)


# ---------------------------------------------
## GB GM

headers = [1,3,5]

from multiprocessing import Process, Manager, TimeoutError, Pool
import numpy as np


# Set the time limit for processing a molecule
time_limit = 15  # in seconds

# Define a function for processing a molecule within the time limit
def process_molecule(smiles, mutants_number, time_limit):
    try:
        start_time = time()
        # Assuming Chem.MolFromSmiles is from RDKit
        original_molecule = Chem.MolFromSmiles(smiles)
        num_atoms = original_molecule.GetNumAtoms()

        # Generate mutants using the specified properties
        input_data = [{"SMILES": smiles, "Num_Mutants": mutants_number, "Max_Atoms": (num_atoms*1.3), "Average_Size": num_atoms, "Size_Stdev": (num_atoms*0.1913)}]
        mutated_smiles = gb_gm_mutator(input_data)  
        elapsed_time = time() - start_time

        # Check if the elapsed time exceeds the time limit
        if elapsed_time > time_limit:
            print(f"Skipping molecule '{smiles}' as it took more than {time_limit} seconds.")
            return None

        # Return the result if successful
        return mutated_smiles

    except Exception as e:
        print(f"Error processing SMILES: {smiles}")
        print(f"Error details: {e}")
        return None

# Initialize lists to store results and pairs
time_taken = []
percent_validity = []
pairs = []

# Initialize headers for CSV files
headers = [1, 3, 5]

# Iterate over each mutant number
for mutants_number in tqdm(headers, desc="Generating mutants"):
    all_mutants = []
    overall_start_time = time()   
    pairs = [] 

    # Generate mutants for all SMILES in the list
    i = 0
    for smiles in smiles_list:
        try:
            i += 1
            print("Method: GB GM Index:", i)
        
            # Use multiprocessing to enforce the time limit
            with Pool(1) as pool:
                result = pool.apply_async(process_molecule, (smiles, mutants_number, 15))
                mutated_smiles = result.get(timeout=15)

            # Check if mutants are generated
            if mutated_smiles:
                # Append the mutated SMILES to the all_mutants list
                all_mutants.extend(mutated_smiles)
                
                # Append the pair of original SMILES and its mutants to the list
                pairs.append((smiles, mutated_smiles))

        except TimeoutError:
            print(f"Skipping molecule '{smiles}' as it took more than 15 seconds.")
            continue


    # End overall time measurement
    overall_end_time = time()
    overall_time_taken = overall_end_time - overall_start_time
    time_taken.append(overall_time_taken)
    
    # Calculate percentual validity for all mutants generated with the current mutant number
    valid_count = sum(IsCorrectSMILES(mutant) for mutant in all_mutants)
    percent_validity.append(valid_count / (mutants_number * len(smiles_list)) * 100)

    # Open a text file in write mode for mutants
    mutants_filename = f'GB_GM_mutants_{mutants_number}.txt'
    mutants_file_path = os.path.join(results_path, mutants_filename)
    with open(mutants_file_path, 'w') as f:
        # Write each mutant to the file
        for mutant in all_mutants:
            f.write(mutant + '\n')

    # Save pairs to file
    pairs_filename = f"GB_GM_pairs_{mutants_number}.txt"
    pairs_file_path = os.path.join(results_path, pairs_filename)
    with open(pairs_file_path, "w") as file:
        for pair in pairs:
            original_smiles, mutants = pair
            file.write("Original SMILES: " + original_smiles + "\n")
            file.write("Mutants:\n")
            for mutant in mutants:
                if mutant is not None:  # Check if mutant is not None
                    file.write(mutant + "\n")
            file.write("\n")  # Add a blank line between pairs

# Create DataFrames for time taken and percentual validity
time_df = pd.DataFrame({"Time Taken (seconds)": time_taken}, index=headers).transpose()
validity_df = pd.DataFrame({"Percentual Validity": percent_validity}, index=headers).transpose()

# Define the filenames for saving CSV files
time_filename = "GB_GM_time_taken.csv"
validity_filename = "GB_GM_percentual_validity.csv"

# Save time_df to CSV with specified path and filename
time_file_path = os.path.join(results_path, time_filename)
time_df.to_csv(time_file_path, index=False)

# Save validity_df to CSV with specified path and filename
validity_file_path = os.path.join(results_path, validity_filename)
validity_df.to_csv(validity_file_path, index=False)

# Print overall time taken
print("Overall Time Taken (seconds):", overall_time_taken)


# ---------------------------------------------
## GB GA 

# Define headers for CSV files
headers = [1, 3, 5]

# Initialize lists to store results
time_taken = []
percent_validity = []

# Set the time limit for processing a molecule
time_limit = 15  # in seconds

# Define a function for processing a molecule within the time limit
def process_molecule(smiles, mutants_number):
    try:
        start_time = time()

        input_data = [{"SMILES": smiles, "Num_Mutants": mutants_number}]

        # Generate mutants using the specified properties
        mutated_smiles = gb_ga_mutator(input_data)

        elapsed_time = time() - start_time

        # Check if the elapsed time exceeds the time limit
        if elapsed_time > time_limit:
            print(f"Skipping molecule '{smiles}' as it took more than {time_limit} seconds.")
            return None

        # Return the result if successful
        return mutated_smiles

    except Exception as e:
        print(f"Error processing SMILES: {smiles}")
        print(f"Error details: {e}")
        return None
    
# Iterate over each mutant number
for mutants_number in tqdm(headers, desc="Generating mutants"):
    
    all_mutants = []
    pairs = []
    overall_start_time = time()    

    # Generate mutants for all SMILES in the list
    i = 0
    for smiles in smiles_list:
        try:
            i += 1
            print("Method: GB GA Index:", i)
        
            # Use multiprocessing to enforce the time limit
            with multiprocessing.Pool(1) as pool:
                result = pool.apply_async(process_molecule, (smiles, mutants_number))
                mutated_smiles = result.get(timeout=time_limit)

            if mutated_smiles is not None:
                # Filter out None values before adding to the list
                all_mutants.extend(filter(None, mutated_smiles))

                if mutated_smiles:
                    # Append the pair of original SMILES and its mutants to the list
                    pairs.append((smiles, mutated_smiles))

        except multiprocessing.TimeoutError:
            print(f"Skipping molecule '{smiles}' as it took more than {time_limit} seconds.")
            continue

    # End overall time measurement
    overall_end_time = time()
    overall_time_taken = overall_end_time - overall_start_time
    time_taken.append(overall_time_taken)
    
    # Calculate percentual validity for all mutants generated with the current mutant number
    valid_count = sum(IsCorrectSMILES(mutant) for mutant in all_mutants)
    percent_validity.append(valid_count / (mutants_number * len(smiles_list)) * 100)

    # Open a text file in write mode
    filename = f'GB_GA_mutants_{mutants_number}.txt'
    file_path = os.path.join(results_path, filename)
    with open(file_path, 'w') as f:
        # Write each mutant to the file
        for mutant in all_mutants:
            f.write(mutant + '\n')

    # Open a text file in write mode for pairs
    filename = f'GB_GA_pairs_{mutants_number}.txt'
    pairs_file_path = os.path.join(results_path, filename)  
    with open(pairs_file_path, "w") as file:
        for pair in pairs:
            original_smiles, mutants = pair
            file.write("Original SMILES: " + original_smiles + "\n")
            file.write("Mutants:\n")
            for mutant in mutants:
                if mutant is not None:  # Check if mutant is not None
                    file.write(mutant + "\n")
            file.write("\n")  # Add a blank line between pairs


# Create DataFrames and save as CSVs (transposed) with specified headers
time_df = pd.DataFrame({"Time Taken (seconds)": time_taken}, index=headers).transpose()
validity_df = pd.DataFrame({"Percentual Validity": percent_validity}, index=headers).transpose()

# Define the filenames for saving CSV files
time_filename = "GB_GA_time_taken.csv"
validity_filename = "GB_GA_percentual_validity.csv"

# Save time_df to CSV with specified path and filename
time_file_path = os.path.join(results_path, time_filename)
time_df.to_csv(time_file_path, index=False)

# Save validity_df to CSV with specified path and filename
validity_file_path = os.path.join(results_path, validity_filename)
validity_df.to_csv(validity_file_path, index=False)

# Print overall time taken
print("Overall Time Taken (seconds):", overall_time_taken)

