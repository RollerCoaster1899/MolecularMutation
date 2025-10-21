import random

def inversion_mutator(array, mutation_rate, num_mutants=1):
    """
    Mutates the input array by inverting a subset of genes based on the mutation rate.

    Parameters:
    - array (list): The input array to be mutated.
    - mutation_rate (float): The probability of a gene subset being inverted.
    - num_mutants (int): The number of mutants to be generated. Default is 1.

    Returns:
    - list: A list containing the mutated arrays.
    """
    mutated_arrays = []

    for _ in range(num_mutants):
        mutated_array = array.copy()  # Create a copy to avoid modifying the original array

        # Perform an inversion only if a randomly generated value is less than the mutation rate
        if random.random() < mutation_rate:
            # Choose a subset of genes to invert (at least two elements)
            subset_size = random.randint(2, len(mutated_array))
            subset_start = random.randint(0, len(mutated_array) - subset_size)

            # Invert the values of the selected subset
            subset_values = mutated_array[subset_start:subset_start + subset_size][::-1]
            mutated_array[subset_start:subset_start + subset_size] = subset_values

        mutated_arrays.append(mutated_array)

    return mutated_arrays