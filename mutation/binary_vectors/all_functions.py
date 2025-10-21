import random  

def swap_mutator(array, mutation_rate):
    """
    Mutates the input array by randomly swapping the positions of two elements
    based on the mutation rate.

    Parameters:
    - array (list): The input array to be mutated.
    - mutation_rate (float): The probability of a swap occurring.

    Returns:
    - list: The mutated array.
    """
    mutated_array = array.copy()  # Create a copy to avoid modifying the original array

    # Perform a swap only if a randomly generated value is less than the mutation rate
    if random.random() < mutation_rate:
        # Select two distinct indices randomly
        index1, index2 = random.sample(range(len(mutated_array)), 2)

        # Swap the values at the selected indices
        mutated_array[index1], mutated_array[index2] = mutated_array[index2], mutated_array[index1]

    return mutated_array

def scramble_mutator(array, mutation_rate):
    """
    Mutates the input array by randomly shuffling a subset of genes based on the mutation rate.

    Parameters:
    - array (list): The input array to be mutated.
    - mutation_rate (float): The probability of a gene subset being shuffled.

    Returns:
    - list: The mutated array.
    """
    mutated_array = array.copy()  # Create a copy to avoid modifying the original array

    # Perform a scramble only if a randomly generated value is less than the mutation rate
    if random.random() < mutation_rate:
        # Choose a subset of genes to shuffle (at least two elements)
        subset_size = random.randint(2, len(mutated_array))
        subset_indices = random.sample(range(len(mutated_array)), subset_size)

        # Shuffle the values of the selected subset
        random.shuffle(subset_indices)
        subset_values = [mutated_array[i] for i in subset_indices]
        for i, index in enumerate(subset_indices):
            mutated_array[index] = subset_values[i]

    return mutated_array

def bitflip_mutator(array, mutation_rate):
    """
    Mutates the input array by flipping bits based on the mutation rate.

    Parameters:
    - array (list): The input array to be mutated.
    - mutation_rate (float): The probability of each bit being flipped.

    Returns:
    - list: The mutated array.
    """
    mutated_array = array.copy()  # Create a copy to avoid modifying the original array

    for i in range(len(mutated_array)):
        if random.random() < mutation_rate:
            mutated_array[i] ^= 1  # Flip the bit using XOR

    return mutated_array

def random_resetting_mutator(array, mutation_rate):
    """
    Mutates the input array by randomly resetting elements from 1 to 0
    based on the mutation rate.

    Parameters:
    - array (list): The input array to be mutated.
    - mutation_rate (float): The probability of each 1 being reset to 0.

    Returns:
    - list: The mutated array.
    """
    mutated_array = array.copy()  # Create a copy to avoid modifying the original array

    for i in range(len(mutated_array)):
        if mutated_array[i] == 1 and random.random() < mutation_rate:
            mutated_array[i] = 0

    return mutated_array


def inversion_mutator(array, mutation_rate):
    """
    Mutates the input array by inverting a subset of genes based on the mutation rate.

    Parameters:
    - array (list): The input array to be mutated.
    - mutation_rate (float): The probability of a gene subset being inverted.

    Returns:
    - list: The mutated array.
    """
    mutated_array = array.copy()  # Create a copy to avoid modifying the original array

    # Perform an inversion only if a randomly generated value is less than the mutation rate
    if random.random() < mutation_rate:
        # Choose a subset of genes to invert (at least two elements)
        subset_size = random.randint(2, len(mutated_array))
        subset_start = random.randint(0, len(mutated_array) - subset_size)

        # Invert the values of the selected subset
        subset_values = mutated_array[subset_start:subset_start + subset_size][::-1]
        mutated_array[subset_start:subset_start + subset_size] = subset_values

    return mutated_array