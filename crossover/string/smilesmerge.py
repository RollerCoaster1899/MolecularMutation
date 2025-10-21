import subprocess
import os
import sys
import shutil

# Define a variable to store the last run number
last_run_number = None

BASE_PATH = "/home/raul-acosta/GitHub/SMILESMerge"
SOURCE_COMPOUNDS_PATH = os.path.join(BASE_PATH, "source_compounds", "Compounds.smi")
MERGED_COMPOUNDS_PATH = os.path.join(BASE_PATH, "merged_compounds")
RUN_SMILES_MERGE_PATH = os.path.join(BASE_PATH, "RunSMILESMerge.py")

def rewrite(smiles_list):
    with open(SOURCE_COMPOUNDS_PATH, 'w') as f:
        for smiles in smiles_list:
            f.write(smiles + "    \n")

def clean_up_merged_folders():
    for folder in os.listdir(MERGED_COMPOUNDS_PATH):
        if folder.startswith("Run_"):
            folder_path = os.path.join(MERGED_COMPOUNDS_PATH, folder)
            shutil.rmtree(folder_path)

def run(number_of_crossovers=1000):
    try:
        subprocess.run([sys.executable, RUN_SMILES_MERGE_PATH, "--number_of_crossovers", str(number_of_crossovers)], stdout=subprocess.PIPE, text=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running subprocess: {e}")
        return None

def extract_smiles_from_latest_run():
    run_folders = [f for f in os.listdir(MERGED_COMPOUNDS_PATH) if f.startswith("Run_")]
    sorted_folders = sorted(run_folders, key=lambda x: int(x.split('_')[1]))

    if sorted_folders:
        latest_run_folder = sorted_folders[-1]
        global last_run_number
        last_run_number = int(latest_run_folder.split("_")[1])

        new_smiles_file = os.path.join(MERGED_COMPOUNDS_PATH, latest_run_folder, "New_SMILES.smi")

        with open(new_smiles_file, 'r') as f:
            merged_smiles = f.read()

        smiles_list = [line.split('\t')[0] for line in merged_smiles.strip().split('\n')]
        return smiles_list
    else:
        return None

def crossover(smiles_list, number_of_crossovers):
    clean_up_merged_folders()  # Clean up before each run
    rewrite(smiles_list)
    run(number_of_crossovers)
    return extract_smiles_from_latest_run()


def smiles_merge_crossover(input_smiles, number_of_mutants):
# Example usage:
    input_smiles = input_smiles
    number_of_mutants = number_of_mutants
    extracted_smiles = crossover(input_smiles, number_of_mutants)

    if extracted_smiles is not None:
        return extracted_smiles
    else:
        return None
    
# Example usage:
input_smiles_list = ["CCc1ccc(N=[N+]=[N-])cc1", "O=C(C)Oc1ccccc1C(=O)O"]
number_of_crossovers = 1  # Change this value as needed
extracted_smiles = smiles_merge_crossover(input_smiles_list, number_of_crossovers)
print(extracted_smiles)