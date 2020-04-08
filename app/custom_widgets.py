import pandas as pd
import ipywidgets as widgets
import numpy as np
from SEIR.seir import SEIR
from SEIR.visualization import visualize_seir_computation

max_simulation_time = widgets.IntText(
    value=300,
    description='Max time')
max_step = widgets.FloatText(
    value=0.5,
    step=0.1,
    description='Max step')
method = widgets.Text(
    value='DOP853',
    description='Method')

compartmentalization = widgets.Checkbox(
    description='Enable Compartmentalization',
    tooltip='Description',
    icon=''  # (FontAwesome names without the `fa-` prefix)
)

population = widgets.IntText(
    value=5e6,
    description='Population size')
incubation_period = widgets.IntText(
    value=3,
    description='Incubation period')
incubation_period = widgets.IntSlider(
    value=3,
    min=1,
    max=14,
    step=1,
    continuous_update=False,
    orientation='horizontal',
    description='Incubation period'
)
infectious_period = widgets.IntText(
    value=7,
    description='Infectious period')
initial_R0 = widgets.FloatText(
    value=2.5,
    step=0.1,
    description='Initial R0')
hospitalization_probability = widgets.FloatText(
    value=2.5,
    step=0.1,
    description='Hospitalization probability')
hospitalization_duration = widgets.IntText(
    value=20,
    description='Hospitalization duration')
hospitalization_lag_from_onset = widgets.IntText(
    value=7,
    description='Hospitalization lag from onset')
icu_probability = widgets.FloatText(
    value=0.01,
    step=0.01,
    description='ICU probability')
icu_duration = widgets.IntText(
    value=10,
    description='ICU duration')
icu_lag_from_onset = widgets.IntText(
    value=11,
    description='ICU lag from onset')
death_probability = widgets.FloatText(
    value=0.1,
    step=0.1,
    description='Death probability')
death_lag_from_onset = widgets.IntText(
    value=25,
    description='Death lag from onset')

style = {'description_width': '150px'}

probabilities = widgets.Dropdown(
    options=[True, False],
    description='Probabilities (True/False)',
    disabled=False,
    indent=False,
    style=style
)
population_susceptible = widgets.FloatText(
    value=0.8,
    step=0.1,
    description='Population susceptible',
    style=style)
population_exposed = widgets.FloatText(
    value=0.15,
    step=0.01,
    description='Population exposed',
    style=style)
population_infected = widgets.FloatText(
    value=0.05,
    step=0.01,
    description='Population infected',
    disabled=False,
    style=style
)

run_simulation = widgets.Button(
    description='Run Simulation',
    button_style='success',
    icon='play'
)


def on_button_clicked(b):
    model = SEIR(population=population.value,
                 incubation_period=incubation_period.value,
                 infectious_period=infectious_period.value,
                 initial_R0=initial_R0.value,
                 hospitalization_probability=hospitalization_probability.value,
                 hospitalization_duration=hospitalization_duration.value,
                 hospitalization_lag_from_onset=hospitalization_lag_from_onset.value,
                 icu_probability=icu_probability.value,
                 icu_duration=icu_duration.value,
                 icu_lag_from_onset=icu_lag_from_onset.value,
                 death_probability=death_probability.value,
                 death_lag_from_onset=death_lag_from_onset.value)
    model.set_initial_state(population_exposed=population_exposed.value,
                            population_infected=population_infected.value,
                            probabilities=probabilities.value)
    model.simulate(max_simulation_time=max_simulation_time.value,
                   method=method.value,
                   max_step=max_step.value)
    # Evaluate the solution
    time = np.arange(0, max_simulation_time.value, 1, dtype=int)
    results = model.evaluate_solution(time)
    # Visualize the results
    visualize_seir_computation(results,
                               compartments=["All"],
                               show_individual_compartments=False)


run_simulation.on_click(on_button_clicked)
