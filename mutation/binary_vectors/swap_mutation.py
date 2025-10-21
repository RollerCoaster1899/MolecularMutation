import random

def swap_mutator(array, mutation_rate, num_mutants=1):
    """
    Mutates the input array by randomly swapping the positions of two elements
    based on the mutation rate.

    Parameters:
    - array (list): The input array to be mutated.
    - mutation_rate (float): The probability of a swap occurring.
    - num_mutants (int): The number of mutants to be generated. Default is 1.

    Returns:
    - list: A list containing the mutated arrays.
    """
    mutated_arrays = []

    for _ in range(num_mutants):
        mutated_array = array.copy()  # Create a copy to avoid modifying the original array

        # Perform a swap only if a randomly generated value is less than the mutation rate
        if random.random() < mutation_rate:
            # Select two distinct indices randomly
            index1, index2 = random.sample(range(len(mutated_array)), 2)

            # Swap the values at the selected indices
            mutated_array[index1], mutated_array[index2] = mutated_array[index2], mutated_array[index1]

        mutated_arrays.append(mutated_array)

    return mutated_arrays