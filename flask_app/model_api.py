import numpy as np
import json

from SEIR.seir import SEIR

simulation_parameters = {'max_simulation_time':'float',
'method':'string',
'max_step':'float'}
model_parameters={"population":'float',
    "incubation_period":'float',
    "population":'float',
    "incubation_period":'float',
    "infectious_period":'float',
    "initial_R0":'float',
    "hospitalization_probability":'float',
    "hospitalization_duration":'float',
    "hospitalization_lag_from_onset":'float',
    "icu_probability":'float',
    "icu_duration":'float',
    "icu_lag_from_onset":'float',
    "death_probability":'float',
    "death_lag_from_onset":'float'}
initial_state_parameters = {
    "population_exposed":'float',
    "population_infected":'float',
    "probabilities":'string'
}

def slice_dict(dictionary, templ):
    subdict={}
    for k in dictionary:
        if k in templ:
            if templ[k] == 'string':
                subdict[k] = dictionary[k]
            elif templ[k] == 'float':
                subdict[k] = float(dictionary[k])
    return(subdict)

def run_simulation(parameters):
    model = SEIR(**slice_dict(parameters,model_parameters))
    model.set_initial_state(**slice_dict(parameters,initial_state_parameters))
    model.simulate(**slice_dict(parameters,simulation_parameters))
    time = np.arange(0, int(parameters['max_simulation_time']), 1, dtype=int)
    results = model.evaluate_solution(time)
    return(results.to_json())
 
if __name__ == '__main__':
    with open('example.json') as json_file:
      data = json.load(json_file)
      run_simulation(data)
