import random

def bitflip_mutator(array, mutation_rate, num_mutants=1):
    """
    Mutates the input array by flipping bits based on the mutation rate.

    Parameters:
    - array (list): The input array to be mutated.
    - mutation_rate (float): The probability of each bit being flipped.
    - num_mutants (int): The number of mutants to be generated. Default is 1.

    Returns:
    - list: A list containing the mutated arrays.
    """
    mutated_arrays = []

    for _ in range(num_mutants):
        mutated_array = array.copy()  # Create a copy to avoid modifying the original array
        for i in range(len(mutated_array)):
            if random.random() < mutation_rate:
                mutated_array[i] ^= 1  # Flip the bit using XOR
        mutated_arrays.append(mutated_array)

    return mutated_arrays