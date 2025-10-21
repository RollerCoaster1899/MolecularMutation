import random

def multi_point_crossover(parent1, parent2, crossover_probability, num_offspring=1):
    """
    Performs Multi-Point Crossover on two parent arrays.

    Parameters:
    - parent1 (list): The first parent array.
    - parent2 (list): The second parent array.
    - crossover_probability (float): The probability of crossover occurring.
    - num_offspring (int): The number of offspring to generate. Default is 1.

    Returns:
    - list: A list containing the generated offspring arrays.
    """
    offspring_arrays = []

    for _ in range(num_offspring):
        # Perform crossover only if a randomly generated value is less than the crossover probability
        if random.random() < crossover_probability:
            # Determine the number of crossover points based on a scaled probability
            num_points = scale_probability_to_points(crossover_probability, len(parent1))

            # Select multiple crossover points
            crossover_points = sorted(random.sample(range(1, min(len(parent1), len(parent2))), num_points))

            # Initialize variables for alternating segments
            segment_start = 0
            current_parent = parent1
            offspring = []

            # Perform multi-point crossover
            for point in crossover_points:
                offspring.extend(current_parent[segment_start:point])
                segment_start = point
                current_parent = parent1 if current_parent is parent2 else parent2

            # Add the remaining segment
            offspring.extend(current_parent[segment_start:])

        else:
            # If no crossover occurs, choose one of the parents as the offspring
            offspring = parent1 if random.random() < 0.5 else parent2

        offspring_arrays.append(offspring)

    return offspring_arrays

def scale_probability_to_points(probability, array_length):
    """
    Scales a probability to a number of crossover points between 2 and the array length.

    Parameters:
    - probability (float): The probability to scale.
    - array_length (int): The length of the array.

    Returns:
    - int: The scaled number of crossover points.
    """
    min_points = 2
    max_points = array_length
    scaled_points = round(min_points + probability * (max_points - min_points))
    return min(max_points, max(min_points, scaled_points))
