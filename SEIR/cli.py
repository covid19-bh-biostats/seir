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
@click.option('--visualize-compartments', type=bool, default=True)
@click.option('--output_file',
              type=str,
              help='Path to output file.')
def main(config_file, contacts_matrix_file, visualize_compartments,
         output_file):
    """Console script for SEIR."""
    # TODO: Handle somehow the creation of an imports function
    # Setup the model
    if not config_file:
        kwargs, initial_state_kwargs, sim_kwargs, restr_info = parse_config_ini(
            f'{WD}/model_configs/finland')
    else:
        kwargs, initial_state_kwargs, sim_kwargs, restr_info = parse_config_ini(
            config_file)

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
    output_file = output_file or os.path.basename(config_file + ".csv")
    results.to_csv(output_file)

    # Visualize the results
    visualize_seir_computation(
        results,
        compartments=kwargs['compartments'],
        restrictions_info=restr_info,
        show_individual_compartments=visualize_compartments)


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
