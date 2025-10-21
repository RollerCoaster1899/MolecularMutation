import random

def one_point_crossover(parent1, parent2, crossover_probability, crossover_point_probability, num_offspring=1):
    """
    Performs One-Point Crossover on two parent arrays.

    Parameters:
    - parent1 (list): The first parent array.
    - parent2 (list): The second parent array.
    - crossover_probability (float): The probability of crossover occurring.
    - crossover_point_probability (float): The probability of selecting the crossover point.
    - num_offspring (int): The number of offspring to generate. Default is 1.

    Returns:
    - list: A list containing the generated offspring arrays.
    """
    offspring_arrays = []

    for _ in range(num_offspring):
        # Perform crossover only if a randomly generated value is less than the crossover probability
        if random.random() < crossover_probability:
            # Select a random crossover point with the specified probability
            if random.random() < crossover_point_probability:
                crossover_point = random.randint(1, min(len(parent1), len(parent2)) - 1)
            else:
                # If the crossover point is not selected, use the midpoint as the default
                crossover_point = len(parent1) // 2

            # Create offspring by swapping tails of parents at the crossover point
            offspring = parent1[:crossover_point] + parent2[crossover_point:]
        else:
            # If no crossover occurs, choose one of the parents as the offspring
            offspring = parent1 if random.random() < 0.5 else parent2

        offspring_arrays.append(offspring)

    return offspring_arrays
