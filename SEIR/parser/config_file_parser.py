import configparser
import itertools
import re

import numpy as np

def _parse_compartments(compartments_str, all_compartments):
    # Check if the input is a single compartment
    if compartments_str in all_compartments:
        return compartments_str

    # 'all' is interpreted as in all compartments
    if re.match('(?i)(all)', compartments_str):
        return all_compartments
    else:
        # Try to interpret the compartments string as a list
        m = re.match(r'\[(?P<compartments>(.+))\]', compartments_str)
        if not m:
            raise ValueError(
                f"[{section_name}]: Compartment list not enclosed in brackets")

        # Different compartments are comma separated, strip any whitespace
        compartments = [c.strip() for c in m.group('compartments').split(',')]

        # See if all the compartments in the infectivity modifier
        # exist
        unknown_compartments = set(compartments).difference(all_compartments)
        if len(unknown_compartments) != 0:
            raise ValueError(
                f"[{section_name}]: Unknown compartments in the infectivity modifier: {unknown_compartments}"
            )
        return compartments


def _parse_restriction_section(section, section_name, all_compartments):
    day_begins = int(section['day-begins'])
    day_ends = int(section['day-ends'])
    info = {
        'begins': day_begins,
        'ends': day_ends,
        'title': " ".join(section_name.split(' ')[1:])
    }
    inf_modifier_str = section['infectivity modifier']
    try:
        # Try to interpret the infectivity modifier as a single float
        inf_modifier = float(inf_modifier_str)

        def restrictions_function(t):
            if day_begins <= t <= day_ends:
                return inf_modifier
            else:
                return 1.0
    except Exception:
        try:
            # Try to interpret the infectivity modifier as a filepath
            raise Exception()
        except Exception:
            # Try to interpret the infectivity modifier as a list of modifications
            # to the matrix
            modifier_str_list = inf_modifier_str.strip().split('\n')
            modifiers = []
            for row in modifier_str_list:
                parts = [part.strip() for part in row.split(':')]
                if len(parts) != 3:
                    raise ValueError(
                        f"{section_name}: Malformed infectivity modifier")

                # Parse the compartments in this modifier
                compartments_1 = _parse_compartments(parts[0],
                                                     all_compartments)
                compartments_2 = _parse_compartments(parts[1],
                                                     all_compartments)

                # Map the compartments to indices
                comp1_idxs = [
                    all_compartments.index(compartment)
                    for compartment in compartments_1
                ]
                comp2_idxs = [
                    all_compartments.index(compartment)
                    for compartment in compartments_2
                ]

                infectivity_matrix_indexes_to_modify = list(itertools.product(comp1_idxs, comp2_idxs))
                try:
                    modifier = float(parts[2])
                except Exception:
                    raise ValueError(
                        f"{section_name}: Malformed infectivity modifier")
                modifier_matrix = np.ones((len(all_compartments), len(all_compartments)))
                modifier_matrix[
                    np.array(infectivity_matrix_indexes_to_modify)[:,0],
                    np.array(infectivity_matrix_indexes_to_modify)[:,1]
                ] *= modifier
                modifier_matrix[
                    np.array(infectivity_matrix_indexes_to_modify)[:,1],
                    np.array(infectivity_matrix_indexes_to_modify)[:,0]
                ] *= modifier

                def restrictions_function(t):
                    if day_begins <= t <= day_ends:
                        return modifier_matrix
                    else:
                        return 1.0

    return restrictions_function, info


def _parse_restriction_sections(config: configparser.ConfigParser,
                                compartments):
    restriction_sections = [
        sec for sec in config.sections()
        if sec.lower().startswith('restriction')
    ]

    restr_funs = []
    restr_info = []
    for rsec in restriction_sections:
        fun, info = _parse_restriction_section(dict(config.items(rsec)), rsec,
                                               compartments)
        restr_funs.append(fun)
        restr_info.append(info)

    if len(restr_funs) == 0:
        return None, None
    elif len(restr_funs) == 1:
        return restr_funs[0], info
    else:

        def restrictions(t):
            modif = 1.0
            for fun in restr_funs:
                modif = np.multiply(modif, fun(t))
            return modif

        return restrictions, info


def parse_config_ini(config_file):
    config = configparser.ConfigParser()
    config.optionxform = str  # Avoids lower-casing of keys
    config.read(config_file)
    kwargs = {}
    for key, value in config.items("model"):
        if "," in value:
            if key == "compartments":
                kwargs[key] = [v.strip() for v in value.split(",")]
            else:
                kwargs[key] = [float(x) for x in value.split(",")]
        else:
            kwargs[key] = float(value)

    # Set default compartment if none given
    if 'compartments' not in kwargs:
        kwargs['compartments'] = ['All']

    # Parse initial state
    initial_state_kwargs = {}
    for key, value in config.items("initial state"):
        if "," in value:
            initial_state_kwargs[key] = [float(x) for x in value.split(",")]
        else:
            try:
                initial_state_kwargs[key] = float(value)
            except Exception:
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
            except Exception:
                if value.lower() in ['yes', 'yau', 'true']:
                    simulation_kwargs[key] = True
                elif value.lower() in ['no', 'nay', 'false']:
                    simulation_kwargs[key] = False
                else:
                    simulation_kwargs[key] = value

    kwargs[
        'restrictions_function'], restrictions_info = _parse_restriction_sections(
            config, kwargs['compartments'])

    return kwargs, initial_state_kwargs, simulation_kwargs, restrictions_info
