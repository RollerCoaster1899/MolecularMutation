import numpy
from conversion.smiles_to_ECFP4 import smiles_to_ECFP4
from conversion.selfies_to_smiles import selfies_to_smiles
from conversion.smiles_to_selfies import smiles_to_selfies
from conversion.ECFP4_to_smiles import ECFP4_to_smiles

smiles = "CCCCC"

# SMILES to ECFP4
ecfp4_mf_p1 = smiles_to_ECFP4(smiles)
print("Original array: " + str(ecfp4_mf_p1) + " Original array sum: " + str(numpy.sum(ecfp4_mf_p1)))

# SMILES to SELFIES
smiles_string = smiles
smiles_list = [smiles, smiles, smiles]

selfies_single = smiles_to_selfies(smiles_string)
selfies_multiple = smiles_to_selfies(smiles_list)

print("SELFIES representation (Single SMILES):", selfies_single)
print("SELFIES representations (Multiple SMILES):", selfies_multiple)

# SELFIES to SMILES
smiles_single = selfies_to_smiles(selfies_single)
smiles_multiple = selfies_to_smiles(selfies_multiple)

print("SMILES representation (Single SELFIES):", smiles_single)
print("SMILES representations (Multiple SELFIES):", smiles_multiple)

# ECFP4 to SMILES
smiles = ECFP4_to_smiles(ecfp4_mf_p1)
print("Converted smiles:", smiles)
print("Type of smiles:", type(smiles))
