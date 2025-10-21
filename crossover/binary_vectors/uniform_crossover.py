import random

def uniform_crossover(parent1, parent2, crossover_probability, bias=None, num_offspring=1):
    """
    Performs Uniform Crossover on two parent arrays.

    Parameters:
    - parent1 (list): The first parent array.
    - parent2 (list): The second parent array.
    - crossover_probability (float): The probability of crossover occurring for each gene.
    - bias (float or None): The bias towards one parent. If None, no bias is applied.
    - num_offspring (int): The number of offspring to generate. Default is 1.

    Returns:
    - list: A list containing the generated offspring arrays.
    """
    offspring_arrays = []

    for _ in range(num_offspring):
        # Perform crossover only if a randomly generated value is less than the crossover probability
        if random.random() < crossover_probability:
            # Generate offspring by flipping a coin for each gene
            offspring = [
                gene_from_parent(parent1[i], parent2[i], bias) for i in range(min(len(parent1), len(parent2)))
            ]
        else:
            # If no crossover occurs, choose one of the parents as the offspring
            offspring = parent1 if random.random() < 0.5 else parent2

        offspring_arrays.append(offspring)

    return offspring_arrays

def gene_from_parent(gene1, gene2, bias=None):
    """
    Selects a gene from the parents with optional bias.

    Parameters:
    - gene1: The gene from the first parent.
    - gene2: The gene from the second parent.
    - bias (float or None): The bias towards one parent. If None, no bias is applied.

    Returns:
    - The selected gene.
    """
    if bias is None:
        # No bias, randomly choose one of the parents
        return gene1 if random.random() < 0.5 else gene2
    else:
        # Apply bias to choose a gene from one parent more often
        return gene1 if random.random() < bias else gene2
