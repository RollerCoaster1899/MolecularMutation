import random

def shuffle_crossover(parent1, parent2, crossover_probability=0.7, num_offspring=1):
    """
    Perform Shuffle Crossover on two parents.

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
        # Randomly select crossover point
        crossover_point = random.randint(0, len(parent1) - 1)

        # Create offspring by shuffling genetic material from both parents in a circular manner
        offspring = []
        for i in range(len(parent1)):
            index = (crossover_point + i) % len(parent1)
            offspring.append(parent1[index] if random.choice([True, False]) else parent2[index])

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