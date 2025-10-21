import random

def partially_mapped_crossover(parent1, parent2, crossover_probability, num_offspring):
    """
    Partially Mapped Crossover (PMX) for permutation-based crossovers.

    Parameters:
    - parent1: The first parent permutation.
    - parent2: The second parent permutation.
    - crossover_probability: Probability of crossover occurring.
    - num_offspring: Number of offspring to generate.

    Returns:
    - offspring: List of offspring permutations.
    """
    offspring = []

    for _ in range(num_offspring):
        if random.random() < crossover_probability:
            # Perform crossover
            crossover_point1 = random.randint(0, len(parent1) - 1)
            crossover_point2 = random.randint(crossover_point1 + 1, len(parent1))

            child = [None] * len(parent1)

            # Copy segment between crossover points from parent1 to child
            child[crossover_point1:crossover_point2] = parent1[crossover_point1:crossover_point2]

            # Map values from parent2 to the corresponding positions in child
            for i in range(crossover_point1, crossover_point2):
                if parent2[i] not in child:
                    index = parent1.index(parent2[i])
                    while child[index] is not None:
                        index = parent1.index(parent2[index])
                    child[index] = parent2[i]

            # Copy remaining values from parent2 to child
            for i in range(len(child)):
                if child[i] is None:
                    child[i] = parent2[i]

            offspring.append(child)
        else:
            # If no crossover, simply copy one of the parents
            offspring.append(list(parent1))

    return offspring
