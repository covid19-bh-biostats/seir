import configparser


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
