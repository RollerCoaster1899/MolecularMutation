import random

def scramble_mutator(array, mutation_rate, num_mutants=1):
    """
    Mutates the input array by randomly shuffling a subset of genes based on the mutation rate.

    Parameters:
    - array (list): The input array to be mutated.
    - mutation_rate (float): The probability of a gene subset being shuffled.
    - num_mutants (int): The number of mutants to be generated. Default is 1.

    Returns:
    - list: A list containing the mutated arrays.
    """
    mutated_arrays = []

    for _ in range(num_mutants):
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

        mutated_arrays.append(mutated_array)

    return mutated_arrays
