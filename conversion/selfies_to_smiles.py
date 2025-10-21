from selfies import decoder, exceptions

def selfies_to_smiles(selfies_input):
    # If the input is a string, convert it to a list with a single element
    if isinstance(selfies_input, str):
        selfies_input = [selfies_input]

    smiles_list = []

    for selfies in selfies_input:
        try:
            smiles = decoder(selfies)
            smiles_list.append(smiles)
        except exceptions.DecoderError as e:
            print(f"Error decoding SELFIES '{selfies}': {e}")
            smiles_list.append(None)

    return smiles_list
