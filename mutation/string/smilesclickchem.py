import subprocess
import os
import sys
import shutil  # Import shutil for folder removal

BASE_PATH = "/home/raul-acosta/GitHub/SMILESClickChem"
SOURCE_COMPOUNDS_PATH = os.path.join(BASE_PATH, "source_compounds", "Compound.smi")
MUTATED_COMPOUNDS_PATH = os.path.join(BASE_PATH, "mutated_compounds")
RUN_SMILES_CLICK_CHEM_PATH = os.path.join(BASE_PATH, "RunSMILESClickChem.py")

def rewrite(smiles):
    with open(SOURCE_COMPOUNDS_PATH, 'w') as f:
        f.write(smiles + "    ")

def clean_up_mutated_folders():
    for folder in os.listdir(MUTATED_COMPOUNDS_PATH):
        if folder.startswith("Run_"):
            folder_path = os.path.join(MUTATED_COMPOUNDS_PATH, folder)
            shutil.rmtree(folder_path)

def run(number_of_mutants=5):
    try:
        subprocess.run([sys.executable, RUN_SMILES_CLICK_CHEM_PATH, "--number_of_mutants", str(number_of_mutants)], stdout=subprocess.PIPE, text=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running subprocess: {e}")
        return None

def extract_smiles():
    run_folders = [f for f in os.listdir(MUTATED_COMPOUNDS_PATH) if f.startswith("Run_")]
    sorted_folders = sorted(run_folders, key=lambda x: int(x.split('_')[1]))

    if sorted_folders:
        latest_run_folder = sorted_folders[-1]
        new_smiles_file = os.path.join(MUTATED_COMPOUNDS_PATH, latest_run_folder, "New_SMILES.smi")

        if os.path.isfile(new_smiles_file):
            with open(new_smiles_file, 'r') as f:
                lines = f.readlines()

            smiles_list = [line.strip().split('\t')[0] for line in lines if line.strip()]
            return smiles_list

    return None

def run_and_extract(smiles, number_of_mutants=5):
    clean_up_mutated_folders()  # Clean up before each run
    rewrite(smiles)
    run(number_of_mutants)
    return extract_smiles()

def smiles_click_chem_mutator(input_smiles, number_of_mutants):
# Example usage:
    input_smiles = input_smiles
    number_of_mutants = number_of_mutants #
    extracted_smiles = run_and_extract(input_smiles, number_of_mutants)

    if extracted_smiles is not None:
        return extracted_smiles
    else:
        return None