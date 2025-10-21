import random

def random_resetting_mutator(array, mutation_rate, num_mutants=1):
    """
    Mutates the input array by randomly resetting elements from 1 to 0
    based on the mutation rate.

    Parameters:
    - array (list): The input array to be mutated.
    - mutation_rate (float): The probability of each 1 being reset to 0.
    - num_mutants (int): The number of mutants to be generated. Default is 1.

    Returns:
    - list: A list containing the mutated arrays.
    """
    mutated_arrays = []

    for _ in range(num_mutants):
        mutated_array = array.copy()  # Create a copy to avoid modifying the original array

        for i in range(len(mutated_array)):
            if mutated_array[i] == 1 and random.random() < mutation_rate:
                mutated_array[i] = 0

        mutated_arrays.append(mutated_array)

    return mutated_arrays