from selfies import encoder, exceptions

def smiles_to_selfies(smiles_input):
    # If the input is a string, convert it to a list with a single element
    if isinstance(smiles_input, str):
        smiles_input = [smiles_input]

    selfies_list = []

    for smiles in smiles_input:
        try:
            selfies = encoder(smiles)
            selfies_list.append(selfies)
        except exceptions.EncoderError as e:
            print(f"Error encoding SMILES '{smiles}': {e}")
            selfies_list.append(None)

    return selfies_list
