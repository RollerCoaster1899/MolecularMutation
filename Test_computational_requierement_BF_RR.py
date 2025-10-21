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
from selfies import encoder, decoder
from conversion.smiles_to_ECFP4 import smiles_to_ECFP4
from mutation.binary_vectors.bitflip_mutation import bitflip_mutator
from mutation.binary_vectors.resetting_mutation import random_resetting_mutator
from mutation.binary_vectors.swap_mutation import swap_mutator
from mutation.binary_vectors.inversion_mutation import inversion_mutator
from mutation.binary_vectors.scramble_mutation import scramble_mutator
from mutation.encoded.smiles_mutation import mutate_smiles
from mutation.encoded.selfies_mutation import mutate_selfies
from mutation.string.smilesclickchem import smiles_click_chem_mutator
from mutation.graph.gb_gm import gb_gm_mutator
from conversion.ECFP4_to_smiles import ECFP4_to_smiles
from itertools import chain
import multiprocessing
from multiprocessing import TimeoutError, Manager, Process


# Load the dataset
df = pd.read_csv("/home/raul-acosta/GitHub/mutation_and_recombination/testing/diversity_and_complexity/comparison_datasets/fda.csv")
smiles_list = df["smiles"].tolist()
smiles_list = random.sample(smiles_list,250)
print(smiles_list)

headers = [5]
# ---------------------------------------------

from multiprocessing import Process, Manager, TimeoutError
import numpy as np

# Initialize lists to store results
time_taken = []
percent_validity = []

# Define a function to process each molecule
def process_molecule(smiles, mutants_number, time_limit, output_list):
    try:
        ecfp4_mf = smiles_to_ECFP4(smiles)
        print("Smiles:", smiles)

        # Generate mutants using the specified properties
        mutated_ECFP4 = bitflip_mutator(ecfp4_mf, 1/len(ecfp4_mf), mutants_number)

        # Convert each mutated ECFP4 to SMILES separately
        mutated_smiles = [ECFP4_to_smiles(mutated) for mutated in mutated_ECFP4]

        # Filter out None values before adding to the list
        output_list.extend(filter(None, mutated_smiles))
    except Exception as e:
        print(f"Error processing SMILES: {smiles}")
        print(f"Error details: {e}")

# Define a function to process a chunk of molecules
# Define a function to process a chunk of molecules
def process_chunk(chunk_index, chunk, mutants_number, time_limit, output_list):
    results = []
    for idx, smiles in enumerate(chunk):
        result = process_molecule(smiles, mutants_number, time_limit, output_list)
        results.append((chunk_index * chunk_size + idx, result))  # Store index along with result
    return results

# Define chunk size
chunk_size = 10

# Iterate over each mutant number
for mutants_number in tqdm(headers, desc="Generating mutants"):
    # Start overall time measurement
    overall_start_time = time()
    all_mutants = []

    # Split smiles_list into chunks
    chunks = [smiles_list[i:i + chunk_size] for i in range(0, len(smiles_list), chunk_size)]

    # Create a manager to share a list among processes
    with Manager() as manager:
        output_list = manager.list()

        # Create a list to store each process
        processes = []

        # Process each chunk in a separate process
        for idx, chunk in enumerate(chunks):
            process = Process(target=process_chunk, args=(idx, chunk, mutants_number, 15, output_list))
            process.start()
            processes.append(process)

        # Wait for each process to finish or timeout
        for process in processes:
            process.join(timeout=155)
            if process.is_alive():
                process.terminate()
                print("Timeout: Process terminated.")

        # Collect results from the shared list
        for results in output_list:
            all_mutants.extend(results)


    # End overall time measurement
    overall_end_time = time()
    overall_time_taken = overall_end_time - overall_start_time
    time_taken.append(overall_time_taken)

    # Calculate percentual validity for all mutants generated with the current mutant number
    valid_count = sum(IsCorrectSMILES(mutant) for mutant in all_mutants)
    percent_validity.append(valid_count / (mutants_number * len(smiles_list)) * 100)

# Create DataFrames and save as CSVs (transposed) with specified headers
time_df = pd.DataFrame({"Time Taken (seconds)": time_taken}, index=headers).transpose()
validity_df = pd.DataFrame({"Percentual Validity": percent_validity}, index=headers).transpose()

time_df.to_csv("BF_time_taken.csv", index=False)
validity_df.to_csv("BF_percentual_validity.csv", index=False)

# Print overall time taken
print("Overall Time Taken (seconds):", overall_time_taken)


# ---------------------------------------------

# Random Resseting
from multiprocessing import Process, Manager, TimeoutError
import numpy as np

# Initialize lists to store results
time_taken = []
percent_validity = []

# Initialize headers for CSV files
headers = [5]

# Define a function to process each molecule
def process_molecule(smiles, mutants_number, time_limit, output_list):
    try:
        ecfp4_mf = smiles_to_ECFP4(smiles)
        print("Smiles:", smiles)

        # Generate mutants using the specified properties
        mutated_ECFP4 = random_resetting_mutator(ecfp4_mf, 1/len(ecfp4_mf), mutants_number)

        # Convert each mutated ECFP4 to SMILES separately
        mutated_smiles = [ECFP4_to_smiles(mutated) for mutated in mutated_ECFP4]

        # Filter out None values before adding to the list
        output_list.extend(filter(None, mutated_smiles))
    except Exception as e:
        print(f"Error processing SMILES: {smiles}")
        print(f"Error details: {e}")

# Define a function to process each molecule
def process_molecule(smiles, mutants_number, time_limit, output_list):
    try:
        ecfp4_mf = smiles_to_ECFP4(smiles)
        print("Smiles:", smiles)

        # Generate mutants using the specified properties
        mutated_ECFP4 = random_resetting_mutator(ecfp4_mf, 1/len(ecfp4_mf), mutants_number)

        # Convert each mutated ECFP4 to SMILES separately
        mutated_smiles = [ECFP4_to_smiles(mutated) for mutated in mutated_ECFP4]

        # Filter out None values before adding to the list
        output_list.extend(filter(None, mutated_smiles))
    except Exception as e:
        print(f"Error processing SMILES: {smiles}")
        print(f"Error details: {e}")

# Define a function to process a chunk of molecules
def process_chunk(chunk, mutants_number, time_limit, output_list):
    for smiles in chunk:
        process_molecule(smiles, mutants_number, time_limit, output_list)

# Define chunk size
chunk_size = 10

# Iterate over each mutant number
for mutants_number in tqdm(headers, desc="Generating mutants"):
    # Start overall time measurement
    overall_start_time = time()
    all_mutants = []

    # Split smiles_list into chunks
    chunks = [smiles_list[i:i + chunk_size] for i in range(0, len(smiles_list), chunk_size)]

    # Create a manager to share a list among processes
    with Manager() as manager:
        output_list = manager.list()

        # Create a list to store each process
        processes = []

        # Process each chunk in a separate process
        for chunk in chunks:
            process = Process(target=process_chunk, args=(chunk, mutants_number, 15, output_list))
            process.start()
            processes.append(process)

        # Wait for each process to finish or timeout
        chunk_n = 0
        for process in processes:
            process.join(timeout=150)
            if process.is_alive():
                process.terminate()
                print("Timeout: Process terminated.")
            chunk_n = chunk_n + 1
            print("index:",chunk_n*chunk_size)  # Print index of each chunk processed

        # Collect results from the shared list
        all_mutants.extend(output_list)

    # End overall time measurement
    overall_end_time = time()
    overall_time_taken = overall_end_time - overall_start_time
    time_taken.append(overall_time_taken)

    # Calculate percentual validity for all mutants generated with the current mutant number
    valid_count = sum(IsCorrectSMILES(mutant) for mutant in all_mutants)
    percent_validity.append(valid_count / (mutants_number * len(smiles_list)) * 100)

# Create DataFrames and save as CSVs (transposed) with specified headers
time_df = pd.DataFrame({"Time Taken (seconds)": time_taken}, index=headers).transpose()
validity_df = pd.DataFrame({"Percentual Validity": percent_validity}, index=headers).transpose()

time_df.to_csv("RR_time_taken.csv", index=False)
validity_df.to_csv("RR_percentual_validity.csv", index=False)

# Print overall time taken
print("Overall Time Taken (seconds):", overall_time_taken)







































# ---------------------------------------------

## GB GM

# Set the time limit for processing a molecule
time_limit = 15  # in seconds

# Define a function for processing a molecule within the time limit
def process_molecule(smiles, mutants_number):
    try:
        start_time = time()
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
    
# Initialize lists to store results
time_taken = []
percent_validity = []

# Initialize headers for CSV files
headers = [5]


# Initialize lists to store results
time_taken = []
percent_validity = []

# Iterate over each mutant number
for mutants_number in tqdm(headers, desc="Generating mutants"):
    
    all_mutants = []
    overall_start_time = time()    

    # Generate mutants for all SMILES in the list
    i = 0
    for smiles in smiles_list:
        try:
            i = i + 1
            print("Method: GB GM Index: ", i)
        
            # Use multiprocessing to enforce the time limit
            with multiprocessing.Pool(1) as pool:
                result = pool.apply_async(process_molecule, (smiles, mutants_number))
                mutated_smiles = result.get(timeout=time_limit)

            if mutated_smiles is not None:
                # Filter out None values before adding to the list
                all_mutants.extend(filter(None, mutated_smiles))

        except TimeoutError:
            print(f"Skipping molecule '{smiles}' as it took more than {time_limit} seconds.")
            continue

    # End overall time measurement
    overall_end_time = time()
    overall_time_taken = overall_end_time - overall_start_time
    time_taken.append(overall_time_taken)
    
    # Calculate percentual validity for all mutants generated with the current mutant number
    valid_count = sum(IsCorrectSMILES(mutant) for mutant in all_mutants)
    percent_validity.append(valid_count / (mutants_number * len(smiles_list)) * 100)

# Create DataFrames and save as CSVs (transposed) with specified headers
time_df = pd.DataFrame({"Time Taken (seconds)": time_taken}, index=headers).transpose()
validity_df = pd.DataFrame({"Percentual Validity": percent_validity}, index=headers).transpose()

time_df.to_csv("GB_GM_time_taken.csv", index=False)
validity_df.to_csv("GB_GM_percentual_validity.csv", index=False)

# Print overall time taken
print("Overall Time Taken (seconds):", overall_time_taken)



# ---------------------------------------------

## SCC
def process_molecule(smiles, mutants_number, all_mutants, lock):
    try:
        # Generate mutants using the specified properties
        scc_smiles = smiles_click_chem_mutator(smiles, mutants_number)
        
        with lock:
            all_mutants.extend(scc_smiles)

    except Exception as e:
        print(f"Error processing SMILES: {smiles}")
        print(f"Error details: {e}")

# Load the dataset
df = pd.read_csv("/home/raul-acosta/GitHub/mutation_and_recombination/testing/diversity_and_complexity/comparison_datasets/fda.csv")
smiles_list = df["smiles"].tolist()

# Initialize lists to store results
time_taken = []
percent_validity = []

# Initialize headers for CSV files
headers = [5]

# Iterate over each mutant number
for mutants_number in tqdm(headers, desc="Generating mutants"):
    # Start overall time measurement
    overall_start_time = time()
    all_mutants = Manager().list()
    lock = Manager().Lock()

    # Generate mutants for all SMILES in the list using multiprocessing
    processes = []
    i = 0
    for smiles in smiles_list:
        i = i + 1
        print("Method: SCC Index: ", i)
        process = Process(target=process_molecule, args=(smiles, mutants_number, all_mutants, lock))
        processes.append(process)
        process.start()

    # Wait for all processes to finish or timeout (15 seconds)
    for process in processes:
        process.join(timeout=15)

    # Terminate any remaining processes
    for process in processes:
        if process.is_alive():
            process.terminate()

    # End overall time measurement
    overall_end_time = time()
    overall_time_taken = overall_end_time - overall_start_time
    time_taken.append(overall_time_taken)

    # Calculate percentual validity for all mutants generated with the current mutant number
    valid_count = sum(IsCorrectSMILES(mutant) for mutant in all_mutants)
    percent_validity.append(valid_count / (mutants_number * len(smiles_list)) * 100)

# Create DataFrames and save as CSVs (transposed) with specified headers
time_df = pd.DataFrame({"Time Taken (seconds)": time_taken}, index=headers).transpose()
validity_df = pd.DataFrame({"Percentual Validity": percent_validity}, index=headers).transpose()

time_df.to_csv("SCC_time_taken.csv", index=False)
validity_df.to_csv("SCC_percentual_validity.csv", index=False)

# Print overall time taken
print("Overall Time Taken (seconds):", overall_time_taken)





## SMILES 

# Initialize lists to store results
time_taken = []
percent_validity = []

# Initialize headers for CSV files
headers = [5]

# Iterate over each mutant number
for mutants_number in tqdm(headers, desc="Generating mutants"):
    # Start overall time measurement
    overall_start_time = time()
    all_mutants = []

    # Generate mutants for all SMILES in the list
    i = 0

    for smiles in smiles_list:
        try:
            i = i + 1
            print("Method: SMILES Index: ", i)
            probability = 1 / len(smiles)

            # Generate mutants using the specified properties
            mutated_smiles = mutate_smiles(smiles, mutants_number, probability, smiles_symbols=None)
            all_mutants.extend(mutated_smiles)

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

# Create DataFrames and save as CSVs (transposed) with specified headers
time_df = pd.DataFrame({"Time Taken (seconds)": time_taken}, index=headers).transpose()
validity_df = pd.DataFrame({"Percentual Validity": percent_validity}, index=headers).transpose()

time_df.to_csv("smiles_time_taken.csv", index=False)
validity_df.to_csv("smiles_percentual_validity.csv", index=False)

# Print overall time taken
print("Overall Time Taken (seconds):", overall_time_taken)


# ---------------------------------------------

# selfies

# Initialize lists to store results
time_taken = []
percent_validity = []

# Initialize headers for CSV files
headers = [5]

# Define a function to process each molecule
def process_molecule(smiles, mutants_number, time_limit, output_list):
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
    except Exception as e:
        print(f"Error processing SMILES: {smiles}")
        print(f"Error details: {e}")

# Iterate over each mutant number
for mutants_number in tqdm(headers, desc="Generating mutants"):
    # Start overall time measurement
    overall_start_time = time()
    all_mutants = []

    # Create a manager to share a list among processes
    with Manager() as manager:
        output_list = manager.list()

        # Create a list to store each process
        processes = []

        i = 0
        # Process each molecule in a separate process
        for smiles in smiles_list:
            i = i + 1
            print("Method: SELFIES Index: ", i)
            
            process = Process(target=process_molecule, args=(smiles, mutants_number, 15, output_list))
            process.start()
            processes.append(process)

        # Wait for each process to finish or timeout
        for process in processes:
            process.join(timeout=15)
            if process.is_alive():
                process.terminate()
                print("Timeout: Process terminated.")

        # Collect results from the shared list
        all_mutants.extend(output_list)

    # End overall time measurement
    overall_end_time = time()
    overall_time_taken = overall_end_time - overall_start_time
    time_taken.append(overall_time_taken)

    # Calculate percentual validity for all mutants generated with the current mutant number
    valid_count = sum(IsCorrectSMILES(mutant) for mutant in all_mutants)
    percent_validity.append(valid_count / (mutants_number * len(smiles_list)) * 100)

# Create DataFrames and save as CSVs (transposed) with specified headers
time_df = pd.DataFrame({"Time Taken (seconds)": time_taken}, index=headers).transpose()
validity_df = pd.DataFrame({"Percentual Validity": percent_validity}, index=headers).transpose()

time_df.to_csv("selfies_time_taken.csv", index=False)
validity_df.to_csv("selfies_percentual_validity.csv", index=False)

# Print overall time taken
print("Overall Time Taken (seconds):", overall_time_taken)


