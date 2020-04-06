import datetime

import numpy as np

from seir.model import SEIR
from seir.visualization import visualize_seir_computation

def main():
    # TODO: Create config-file parsing
    # TODO: Handle reading contacts matrix from a file if provided
    # TODO: Handle somehow the creation of a restrictions function
    # TODO: Handle somehow the creation of an imports function

    # Setup the model
    model = SEIR(incubation_period=3,
                 infectious_period=7,
                 initial_R0=2.3,
                 hospitalization_probability=[0.01, 0, 0, 0, 0.2],
                 hospitalization_duration=21,
                 hospitalization_lag_from_onset=6,
                 icu_probability=0.001,
                 icu_duration=7,
                 icu_lag_from_onset=21,
                 death_probability=0.1,
                 death_lag_from_onset=27,
                 compartments=['G1', 'G2', 'G3', 'G4', 'G5'],
                 population=[2.5e6, 1e6, 3e5, 5e4, 4e5])

    # Setup initial state
    model.set_initial_state(population_susceptible=0.99,
                            population_exposed=0.005,
                            population_infected=0.005,
                            probabilities=True)

    # Simulate up to 200 days
    model.simulate(200)

    # Evaluate the solution
    time = np.arange(0, 200, 1, dtype=int)
    results = model.evaluate_solution(time)

    # Save data
    results.to_csv('outfile.csv')

    # Visualize the results
    visualize_seir_computation(results,
                               show_individual_compartments = True)


if __name__ == '__main__':
    main()
