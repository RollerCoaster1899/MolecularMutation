def whole_arithmetic_recombination(parent1, parent2, alpha, num_offspring=2):
    """
    Performs Whole Arithmetic Recombination on two parent arrays.

    Parameters:
    - parent1 (list): The first parent array.
    - parent2 (list): The second parent array.
    - alpha (float): The weighting factor for the recombination.
    - num_offspring (int): The number of offspring to generate. Default is 2.

    Returns:
    - list: A list containing the generated offspring arrays.
    """
    offspring_arrays = []

    for _ in range(num_offspring):
        child = [(alpha * gene1 + (1 - alpha) * gene2) for gene1, gene2 in zip(parent1, parent2)]
        offspring_arrays.append(child)

    return offspring_arrays