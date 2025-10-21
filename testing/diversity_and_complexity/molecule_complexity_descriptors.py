from rdkit import Chem
from rdkit.Chem import Descriptors

def calculate_molecular_weight(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")
    
    molecular_weight = Descriptors.MolWt(mol)
    return molecular_weight

def calculate_fraction_sp3(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")

    num_sp3_carbons = Descriptors.FractionCSP3(mol)
    return num_sp3_carbons

def calculate_fraction_chiral_carbons(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")

    num_chiral_carbons = Descriptors.FractionCSP3(mol)
    return num_chiral_carbons

def calculate_number_of_rings(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")

    num_rings = Descriptors.RingCount(mol)
    return num_rings

def calculate_number_of_heteroatoms(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")

    num_heteroatoms = Descriptors.NumHeteroatoms(mol)
    return num_heteroatoms



