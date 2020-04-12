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


def test_restrictions_content(finland_config_with_restrictions):
    restrictions = finland_config_with_restrictions.restrictions
    assert len(restrictions) == 2
    assert restrictions[0]['begins'] == 20
    assert restrictions[0]['ends'] == 180
    assert restrictions[0]['title'] == 'everyone-begin-careful'

    assert restrictions[1]['begins'] == 40
    assert restrictions[1]['ends'] == 100
    assert restrictions[1]['title'] == 'school-closure'
