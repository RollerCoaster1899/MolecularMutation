import random

def order_based_crossover(parent1, parent2, crossover_probability, num_offspring):
    """
    Perform Order-based Crossover (OX2) on two parents.

    Parameters:
    - parent1: List representing the first parent
    - parent2: List representing the second parent
    - crossover_probability: Probability of performing crossover (default is 0.7)
    - num_offspring: Number of offspring to generate (default is 1)

    Returns:
    - offspring: List of lists representing the new offspring
    """

    # Ensure the input lists have the same length
    if len(parent1) != len(parent2):
        raise ValueError("Parents must have the same length")

    # Helper function for crossover operation
    def perform_crossover(parent1, parent2):
        point1, point2 = sorted(random.sample(range(len(parent1)), 2))
        offspring = [-1] * len(parent1)
        offspring[point1:point2 + 1] = parent1[point1:point2 + 1]
        pointer_offspring = pointer_parent2 = point2 + 1
        pointer_parent1 = (point2 + 1) % len(parent1)

        while offspring.count(-1) > 0:
            if offspring[pointer_offspring % len(parent1)] == -1:
                if parent2[pointer_parent2 % len(parent2)] not in offspring:
                    offspring[pointer_offspring % len(parent1)] = parent2[pointer_parent2 % len(parent2)]
                    pointer_offspring += 1
                pointer_parent2 += 1
            else:
                pointer_offspring += 1

        return offspring

    # Perform crossover for each offspring
    offspring_list = []
    for _ in range(num_offspring):
        if random.random() < crossover_probability:
            offspring_list.append(perform_crossover(parent1, parent2))
        else:
            # If crossover doesn't occur, choose one of the parents as the offspring
            offspring_list.append(random.choice([parent1, parent2]))

    return offspring_list

