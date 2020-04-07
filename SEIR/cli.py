"""Console script for SEIR."""
import os
import sys
import click
import numpy as np

from SEIR.parser.config_file_parser import parse_config_ini
from SEIR.seir import SEIR
from SEIR.visualization import visualize_seir_computation

WD = os.path.dirname(__file__)


@click.command()
@click.option('--config_file',
              type=click.Path(exists=True),
              help='Path to config ini file.')
@click.option('--contacts_matrix_file',
              type=click.Path(exists=True),
              help='Path to contact matrix file')
@click.option('--output_file',
              type=str,
              help='Path to output file.')
def main(config_file, contacts_matrix_file, output_file):
    """Console script for SEIR."""
    # TODO: Handle somehow the creation of a restrictions function
    # TODO: Handle somehow the creation of an imports function
    # Setup the model
    if not config_file:
        kwargs = parse_config_ini(f"{WD}/model_configs/sample.ini")
    else:
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
    results.to_csv(output_file)

    # Visualize the results
    visualize_seir_computation(results,
                               compartments=kwargs['compartments'],
                               show_individual_compartments=True)


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
