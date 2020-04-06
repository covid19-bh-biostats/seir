import datetime
import argparse
import sys

import numpy as np

from seir.model import SEIR
from seir.visualization import visualize_seir_computation


def main(contacts_matrix_file, output_file):
    # TODO: Create config-file parsing
    # TODO: Handle somehow the creation of a restrictions function
    # TODO: Handle somehow the creation of an imports function
    # Setup the model
    kwargs = dict(
        incubation_period=3,
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
        population=[2.5e6, 1e6, 3e5, 5e4, 4e5],
    )
    if contacts_matrix_file:
        with open(contacts_matrix_file) as contacts_matrix_file:
            kwargs["contacts_matrix"] = np.fromfile(contacts_matrix_file, dtype=np.int32)

    model = SEIR(**kwargs)

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
    results.to_csv(output_file)

    # Visualize the results
    visualize_seir_computation(results,
                               show_individual_compartments = True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Modeling epidemics using the SEIR model')
    parser.add_argument('-c', dest="contacts_matrix_file", type=str, help='File path to contact matrix file')
    parser.add_argument('-o', dest="output_file", default="outfile.csv", type=str, help='Output file name')
    args = parser.parse_args(sys.argv[1:])
    main(args.contacts_matrix_file, args.output_file)
