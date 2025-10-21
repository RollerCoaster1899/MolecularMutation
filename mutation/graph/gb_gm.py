import csv

def write_to_csv(output_file_path, data):
    # Specify the fieldnames for the CSV file
    fieldnames = ["SMILES", "Num_Mutants", "Max_Atoms", "Average_Size", "Size_Stdev"]

    # Write data to the CSV file
    with open(output_file_path, "w", newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        # Write header
        writer.writeheader()
        # Write data rows
        for row in data:
            writer.writerow(row)


import subprocess
import os 
def run(script_path):
    try:
        # Get the directory of the script
        script_directory = os.path.dirname(os.path.abspath(script_path))

        # Set the working directory to the script directory
        os.chdir(script_directory)

        # Run the script as a subprocess
        subprocess.run(["python3", script_path], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error: {e}")
    finally:
        # Restore the working directory to the original
        os.chdir(os.path.dirname(os.path.abspath(__file__)))

def extract_smiles_from_file(file_path):
    """
    Extracts the list of SMILES from a saved file.
    
    Parameters:
    - file_path (str): The path to the file containing SMILES.

    Returns:
    - list: A list of SMILES.
    """
    smiles_list = []
    
    try:
        with open(file_path, "r") as file:
            smiles = ""
            for line in file:
                line = line.strip()
                # Concatenate lines until a complete SMILES is formed
                if line.endswith("\\"):
                    smiles += line[:-1]
                else:
                    smiles += line
                    smiles_list.append(smiles)
                    smiles = ""
    except FileNotFoundError:
        print(f"File not found: {file_path}")
    
    return smiles_list

def gb_gm_mutator(input_data):
    script_path = "/home/raul-acosta/GitHub/GB_BM/GB_GM_mutator.py"
    output_file_path = "/home/raul-acosta/GitHub/GB_BM/Compound.csv"
    output_mutator_path = "/home/raul-acosta/GitHub/GB_BM/mutated_compounds/mutated_SMILES.txt"

    write_to_csv(output_file_path, input_data)
    run(script_path)
    mutated_compounds = extract_smiles_from_file(output_mutator_path)

    return mutated_compounds

