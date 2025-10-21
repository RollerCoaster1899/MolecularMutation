import subprocess
import numpy as np
import os
import torch
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def generate_fingerprint(smiles, max_bit=2048, radius=2):
    try:
        mol = Chem.MolFromSmiles(smiles)
        fingerprint = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius, nBits=max_bit, useFeatures=True)

        fingerprint_array = np.zeros(max_bit, dtype=int)
        for pos in fingerprint.GetOnBits():
            fingerprint_array[pos] = 1

        return fingerprint_array
    except:
        return None
    
def calculate_positions_from_fingerprint(fingerprint):
    return [pos for pos in range(len(fingerprint)) if fingerprint[pos] == 1]

def calculate_smiles_from_positions(positions, predict_script_path):
    try:
        positions_str = ' '.join(map(str, positions))
        command = f"python3 {predict_script_path} --fp='ECFP4' --model_type='smiles' --input='{positions_str}'"
        output = subprocess.check_output(command, shell=True, universal_newlines=True)

        result_line = [line for line in output.split('\n') if line.startswith("Result:")][0]
        result_smiles = result_line.split("Result: ")[1].replace(" ", "").strip()

        return result_smiles
    except Exception as e:
        print(f"Error occurred: {str(e)}")
        return None

def ECFP4_to_smiles(ecfp4_array):
    # Calculate positions from ECFP4 fingerprint
    positions = calculate_positions_from_fingerprint(ecfp4_array)

    # Save the current working directory
    original_working_directory = os.getcwd()

    # Change the working directory to the desired path
    new_working_directory = '/home/raul-acosta/GitHub/MolForge'
    os.chdir(new_working_directory)

    # Specify the path to your predict.py script
    predict_script_path = os.path.join('/home/raul-acosta/GitHub/MolForge', 'predict.py')

    try:
        # Calculate SMILES from positions
        predicted_smiles = calculate_smiles_from_positions(positions, predict_script_path)
    finally:
        # Clear CUDA cache after subprocess call
        torch.cuda.empty_cache()

    # Return to the original working directory
    os.chdir(original_working_directory)

    return predicted_smiles if predicted_smiles is not None else None
