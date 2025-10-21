import random

def davis_order_crossover(parent1, parent2):
    """
    Perform Davis' Order Crossover (OX1) on two parents.

    Parameters:
    - parent1: List representing the first parent (permutation)
    - parent2: List representing the second parent (permutation)

    Returns:
    - offspring1: List representing the first offspring
    - offspring2: List representing the second offspring
    """

    # Ensure the input lists have the same length
    if len(parent1) != len(parent2):
        raise ValueError("Parents must have the same length")

    # Choose two random crossover points
    point1, point2 = sorted(random.sample(range(len(parent1)), 2))

    # Create the first offspring by copying the segment between crossover points from the first parent
    offspring1 = parent1[point1:point2]

    # Initialize pointers for the remaining positions in the second parent
    pointer_parent2 = point2 % len(parent2)
    pointer_offspring1 = point2 % len(parent1)

    # Fill in the remaining positions in the first offspring using genetic material from the second parent
    while len(offspring1) < len(parent1):
        if parent2[pointer_parent2] not in offspring1:
            offspring1.insert(pointer_offspring1 % len(parent1), parent2[pointer_parent2])
            pointer_offspring1 += 1
        pointer_parent2 = (pointer_parent2 + 1) % len(parent2)

    # Repeat the process for the second offspring with the parents' roles reversed
    offspring2 = davis_order_crossover(parent2, parent1)

    return offspring1, offspring2

