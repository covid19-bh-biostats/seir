import os
from dataclasses import dataclass

import pytest

from SEIR.parser.config_file_parser import parse_config_ini

WD = os.path.dirname(__file__)


@dataclass
class FullConfig:
    model: dict
    initial_state: dict
    simulation: dict
    restrictions: dict


@pytest.fixture(scope='function')
def finland_config_with_restrictions(rootdir):
    """
    parses test_files_finland_with_restrictions and provides a FullConfig object
    """
    model, initial_state, simulation, restrictions = parse_config_ini(f'{rootdir}/test_files/finland_with_restrictions')

    return FullConfig(model, initial_state, simulation, restrictions)


def test_matching_lengths(finland_config_with_restrictions):
    """
    Verifies that indeed everything in the model config file was parsed
    """
    assert len(finland_config_with_restrictions.model) == 14
    assert len(finland_config_with_restrictions.initial_state) == 3
    assert len(finland_config_with_restrictions.simulation) == 3
    assert len(finland_config_with_restrictions.restrictions) == 2


def test_model_types(finland_config_with_restrictions):
    model = finland_config_with_restrictions.model
    assert isinstance(model['compartments'], list)
    assert all(isinstance(comp, str) for comp in model['compartments'])
    assert isinstance(model['population'], list)
    assert all(isinstance(comp, float) for comp in model['population'])
    assert isinstance(model['incubation_period'], float)
    assert isinstance(model['infectious_period'], float)
    assert isinstance(model['initial_R0'], float)
    assert isinstance(model['hospitalization_probability'], list)
    assert isinstance(model['hospitalization_duration'], float)
    assert isinstance(model['hospitalization_lag_from_onset'], float)
    assert isinstance(model['death_probability'], list)
    assert isinstance(model['death_lag_from_onset'], float)


def test_initial_state_types(finland_config_with_restrictions):
    initial_state = finland_config_with_restrictions.initial_state
    assert isinstance(initial_state['probabilities'], bool)
    assert isinstance(initial_state['population_exposed'], list)
    assert isinstance(initial_state['population_infected'], list)


def test_simulation_types(finland_config_with_restrictions):
    simulation = finland_config_with_restrictions.simulation
    assert isinstance(simulation['max_simulation_time'], float)
    assert isinstance(simulation['method'], str)
    assert isinstance(simulation['max_step'], float)


def test_restrictions_types(finland_config_with_restrictions):
    restrictions = finland_config_with_restrictions.restrictions
    assert isinstance(restrictions, list)
    assert isinstance(restrictions[0], dict)
    assert isinstance(restrictions[1], dict)
