import argparse
import configparser

import numpy as np

from seir.model import SEIR
from seir.visualization import visualize_seir_computation


def parse_config_ini(config_file):
    config = configparser.ConfigParser()
    config.optionxform = str  # Avoids lower-casing of keys
    config.read(config_file)
    kwargs = {}
    for key, value in config.items("model"):
        if "," in value:
            if key == "compartments":
                kwargs[key] = value.split(",")
            else:
                kwargs[key] = [float(x) for x in value.split(",")]
        else:
            kwargs[key] = float(value)
    return kwargs


def main(config_file, contacts_matrix_file, output_file):
    # TODO: Handle somehow the creation of a restrictions function
    # TODO: Handle somehow the creation of an imports function
    # Setup the model
    kwargs = parse_config_ini(config_file)

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
    output_file = output_file or config_file + ".csv"
    results.to_csv(output_file)

    # Visualize the results
    visualize_seir_computation(results,
                               show_individual_compartments=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Modeling epidemics using the SEIR model')
    parser.add_argument('config_file', type=str, help='File path to config ini file')
    parser.add_argument('-c', dest="contacts_matrix_file", type=str, help='File path to contact matrix file')
    parser.add_argument('-o', dest="output_file", type=str, help='Output file name')
    args = parser.parse_args()
    main(args.config_file, args.contacts_matrix_file, args.output_file)
