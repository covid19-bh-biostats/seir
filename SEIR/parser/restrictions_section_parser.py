import configparser
import itertools
import re
import pathlib
from typing import Callable, Dict, List, Optional, Text, Tuple, Union

import numpy as np


def _parse_compartments(compartments_str: Text,
                        all_compartments: List[Text]) -> List[Text]:
    """
    Converts a string representing one or multiple compartments
    to a list of compartments. 'all' is converted to list of
    all compartment names.

    Input
    -----
    compartments_str: Text
        A string representing the compartments. If a list,
        should be comma-separated and enclosed in brackets, e.g.,
        - "0-4", or
        - "[ 0-4, 5-10, 60-64 ], or
        - "all"
    all_compartments: List[Text]
        A list of the names of all compartments in the model

    Returns
    -------
    List[Text]
        compartments_str converted to a list of compartment names

    Raises
    ------
    ValueError
        If the input cannot be interpreted as a list of compartments.
    """
    # Check if the input is a single compartment
    if compartments_str in all_compartments:
        return [compartments_str]

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
            raise ValueError((f"[{section_name}]: Unknown compartments in the "
                              f"infectivity modifier: {unknown_compartments}"))
        return compartments


def _get_infectivity_modifier_from_file(path: Text, num_all_compartments: int
                                        ) -> np.ndarray:
    """
    Loads infecticity modifier (matrix) from a CSV file.

    Input
    -----
    path: Text
        Path to the infectivity modifier file, CSV file. Should
        contain a matrix of shape
        (num_all_compartments, num_all_compartments).
    num_all_compartments: int
        Number of compartments in the model

    Returns
    -------
    np.ndarray
        The modifier matrix of the infectivity rate

    Raises
    ------
    ValueError
        If the file cannot be found or the matrix has wrong shape
    """
    if pathlib.Path(path).is_file():
        modifier_matrix = np.loadtxt(path)

        # Check that the shape is correct
        if modifier_matrix.shape !=\
           (num_all_compartments, num_all_compartments):
            raise ValueError((f"[{section_name}]: Infectivity modifier "
                              f"matrix in {file} has shape: "
                              f"{modifier_matrix.shape}, expected "
                              f"({len(all_compartments)}, "
                              f"{len(all_compartments)})"))

    else:
        raise ValueError((f"[{section_name}]: Infectivity modifier file "
                          f"{path} does not exist."))

    return modifier_matrix


def _parse_infectivity_modifier_matrix_definition_single_line(
        line: Text,
        all_compartments: List[Text]) -> Tuple[List[int], List[int], float]:
    """
    Parses a single line defining modifications to the infectivity matrix.

    Input
    -----
    line: Text
        Line defining the modifications to the infectivity matrix. Should be
        a colon-separated list (without brackets) where two first elements
        define the compartments this modification relates to and the least
        element should be a float defining the prefactor for the element
        of the infectivity matrix
    all_compartments: List[Text]
        List of all defined compartments

    Returns
    -------
    i: List[int]
        Column indices for the modifier matrix.
    j: List[int]
        Row indices for the modifier matrix.
    modifier: float
        The value with which to multiple the corresponding
        elements of the modifier matrix

    Raises
    ------
    ValueError
        Raised if the modifier line is invalid
    """
    # Check that we have three parts, separated by ':'
    parts = [part.strip() for part in line.split(':')]
    if len(parts) != 3:
        raise ValueError(f"{section_name}: Malformed infectivity modifier")

    # Parse the compartments in this modifier
    compartments_1 = _parse_compartments(parts[0], all_compartments)
    compartments_2 = _parse_compartments(parts[1], all_compartments)

    # Map the compartments to indices
    comp1_idxs = [
        all_compartments.index(compartment) for compartment in compartments_1
    ]
    comp2_idxs = [
        all_compartments.index(compartment) for compartment in compartments_2
    ]

    infectivity_matrix_indexes_to_modify = np.array(
        list(itertools.product(comp1_idxs, comp2_idxs)) +
        list(itertools.product(comp2_idxs, comp1_idxs)))

    try:
        modifier = float(parts[2])
    except Exception:
        raise ValueError((f"{section_name}: Malformed infectivity modifier "
                          "{parts[2]}, expected a float"))
    return infectivity_matrix_indexes_to_modify[:, 0],\
        infectivity_matrix_indexes_to_modify[:, 1], modifier


def _parse_infectivity_modifier_matrix_definition_string(
        infectivity_modifier_string: Text,
        all_compartments: List[Text]) -> np.ndarray:
    """
    Parses the infectivity rate/restriction matrix definition string

    Input
    -----
    infectivity_modifier_string: Text
    all_compartments: List[Text]

    Returns
    -------
    np.ndarray: The prefactor matrix for element-wise modification of the
        infectivity rate matrix
    """
    # Generate the initial modifier matrix
    modifier_matrix = np.ones((len(all_compartments), len(all_compartments)))

    # Loop over all the matrix definitions
    modifier_str_list = infectivity_modifier_string.strip().split('\n')
    for row in modifier_str_list:
        idx_i, idx_j, modifier =\
            _parse_infectivity_modifier_matrix_definition_single_line(
                row, all_compartments
            )
        modifier_matrix[idx_i, idx_j] *= modifier
    return modifier_matrix


def _parse_restriction_section(
        section: Dict[Text, Text], section_name: Text,
        all_compartments: List[Text]
) -> Tuple[Callable[[float], Union[float, np.
                                   ndarray]], Dict[Text, Union[Text, int]]]:
    """
    Parses a single [restriction LABEL] section of the config file.

    Input
    -----
    section: Dict[Text, Text]
        Contents of the section in the config file
        as a dictionary.
    section_name: Text
        Name of the section
    all_compartments: List[Text]
        A list of all defined compartment names

    Returns
    -------
    restrictions_function: Callable(float) -> float/np.ndarray
        A function with a time (in days) as its argument. Returns
        the modification factor to the infectivity matrix.
    info: Dict[Text, Union[Text, int]]
        A dictionary with information regarding the restriction:

        begins
            The day the restriction starts
        ends
            The day the restriction ends
        title
            Name of the restriction (from the config file section title)

    Raises
    ------
    ValueError
        If the restriction section cannot be parsed.
    """
    num_all_compartments = len(all_compartments)

    # Setup the info dictionary
    day_begins = int(section['day-begins'])
    day_ends = int(section['day-ends'])
    info = {
        'begins': day_begins,
        'ends': day_ends,
        'title': " ".join(section_name.split(' ')[1:])
    }

    # Parse the infectivity modifier
    inf_modifier_str = section['infectivity modifier']
    inf_modifier = None
    try:
        # Try to interpret the infectivity modifier as a single float
        inf_modifier = float(inf_modifier_str)
    except Exception:
        ...

    if not inf_modifier:  # Is not a float
        # Try to interpret the infectivity modifier as a filepath
        if inf_modifier_str.strip().startswith('file://'):
            path = inf_modifier_str.strip().replace('file://', '')
            inf_modifier = _get_infectivity_modifier_from_file(
                path, num_all_compartments)
        else:
            # Try to interpret the infectivity modifier string as
            # a list of modifications to the matrix
            inf_modifier =\
                _parse_infectivity_modifier_matrix_definition_string(
                    inf_modifier_str, all_compartments)

    def restrictions_function(t: float) -> Union[float, np.ndarray]:
        if day_begins <= t <= day_ends:
            return inf_modifier
        else:
            return 1.0

    return restrictions_function, info


def parse_restriction_sections(
        config: configparser.ConfigParser, compartments: List[Text]
) -> Tuple[Optional[Callable[[float], Union[float, np.ndarray]]], Optional[
        List[Dict[Text, Union[Text, int]]]]]:
    """
    Parses all "[restriction NAME]" sections from the configparser to
    a function and list of infos on the restrictions.

    Input
    -----
    config: configparser.ConfigParser
        The parser
    compartments: List[Text]
        List of all compartments

    Returns
    -------
    restrictions_function(float) -> float/np.ndarray
        A function of time returning the element-wise prefactors
        for the infectivity rate matrix
    restriction_infos: List
        List of infos on the restrictions. The infos are dictionaries
        with the following keys

        begins
            The day the restriction starts
        ends
            The day the restriction ends
        title
            Name of the restriction (from the config file section title)

    Raises
    ------
    ValueError
        If there is an error in the restrictions definition
    """

    # Get the section names for restriction definitions
    restriction_sections = [
        sec for sec in config.sections()
        if sec.lower().startswith('restriction')
    ]

    # Parse all restriction definitions
    restr_funs = []
    restr_info = []
    for rsec in restriction_sections:
        fun, info = _parse_restriction_section(dict(config.items(rsec)), rsec,
                                               compartments)
        restr_funs.append(fun)
        restr_info.append(info)

    # Return the restriction function and infos on the restrictions
    if len(restr_funs) == 0:
        return None, None
    elif len(restr_funs) == 1:
        return restr_funs[0], restr_info[0]
    else:

        def restrictions(t):
            modif = 1.0
            for fun in restr_funs:
                modif = np.multiply(modif, fun(t))
            return modif

        return restrictions, restr_info

