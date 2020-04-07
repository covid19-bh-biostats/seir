import argparse
import configparser
import sys

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

    # Parse initial state
    initial_state_kwargs = {}
    for key, value in config.items("initial state"):
        if "," in value:
            initial_state_kwargs[key] = [float(x) for x in value.split(",")]
        else:
            try:
                initial_state_kwargs[key] = float(value)
            except Exception as e:
                if value.lower() in ['yes', 'yau', 'true']:
                    initial_state_kwargs[key] = True
                elif value.lower() in ['no', 'nay', 'false']:
                    initial_state_kwargs[key] = False
                else:
                    raise ValueError(f"Non-float value for {key}: {value}")

    # Parse simulation arguments
    simulation_kwargs = {}
    for key, value in config.items("simulation"):
        if "," in value:
            simulation_kwargs[key] = [float(x) for x in value.split(",")]
        else:
            try:
                simulation_kwargs[key] = float(value)
            except Exception as e:
                if value.lower() in ['yes', 'yau', 'true']:
                    simulation_kwargs[key] = True
                elif value.lower() in ['no', 'nay', 'false']:
                    simulation_kwargs[key] = False
                else:
                    simulation_kwargs[key] = value
    return kwargs, initial_state_kwargs, simulation_kwargs


def _main_core(config_file, contacts_matrix_file, output_file,
               visualize_compartments):
    # TODO: Handle somehow the creation of a restrictions function
    # TODO: Handle somehow the creation of an imports function
    # Setup the model
    kwargs, initial_state_kwargs, sim_kwargs = parse_config_ini(config_file)

    if contacts_matrix_file:
        with open(contacts_matrix_file) as contacts_matrix_file:
            kwargs["contacts_matrix"] = np.loadtxt(contacts_matrix_file)

    model = SEIR(**kwargs)

    # Setup initial state
    model.set_initial_state(**initial_state_kwargs)

    # Simulate up to 200 days
    model.simulate(**sim_kwargs)

    # Evaluate the solution
    time = np.arange(0, sim_kwargs['max_simulation_time'], 1, dtype=int)
    results = model.evaluate_solution(time)

    # Save data
    results.to_csv(output_file)

    # Visualize the results
    visualize_seir_computation(
        results,
        compartments=kwargs['compartments'],
        show_individual_compartments=visualize_compartments)


def main():
    parser = argparse.ArgumentParser(
        description='Modeling epidemics using the SEIR model')
    parser.add_argument('config_file',
                        type=str,
                        help='File path to config ini file')
    parser.add_argument('-c',
                        dest="contacts_matrix_file",
                        type=str,
                        help='File path to contact matrix file')
    parser.add_argument('--visualize-compartments',
                        action='store_true',
                        help='Visualize dynamics of individual compartments.')
    parser.add_argument('-o',
                        dest="output_file",
                        default="outfile.csv",
                        type=str,
                        help='Output file name')
    args = parser.parse_args()
    _main_core(args.config_file, args.contacts_matrix_file, args.output_file,
               args.visualize_compartments)


if __name__ == '__main__':
    main()
