"""Main module."""
import itertools
from typing import Any, Callable, List, Optional, Union, Text

import numpy as np
from scipy.integrate import solve_ivp, cumtrapz, trapz
import pandas as pd


class SEIR:
    """
    Implementation of a SEIR model
    """

    def __init__(
        self,
        *,
        incubation_period: Union[int, float, np.ndarray],
        infectious_period: Union[int, float, np.ndarray],
        initial_R0: Union[int, float],
        hospitalization_probability: Union[float, np.ndarray],
        hospitalization_duration: Union[float, int],
        hospitalization_lag_from_onset: Union[float, int],
        icu_probability: Union[float, int],
        icu_duration: Union[float, int],
        icu_lag_from_onset: Union[float, np.ndarray],
        death_probability: Union[float, np.ndarray],
        death_lag_from_onset: Union[float, int],
        population: Union[float, np.ndarray],
        compartments: Optional[List[Any]] = None,
        contacts_matrix: Optional[np.ndarray] = None,
        restrictions_function: Optional[Callable[[float], Union[float, np.ndarray]]] = None,
        imported_cases_function: Optional[Callable] = None
    ):
        """
        Initializes the SEIR models parameters and computes the infectivity
        rate from the contacts matrix, R0, and infective_duration. Supports
        compartmentalizing the population to, e.g., age-groups.

        Keyword arguments
        -----------------
        incubation_period: Union[int, float, np.ndarray]
            Incubation period of the disease in days. If an array,
            it should be the incubation period for each population compartment.
        infectious_period: Union[int, float, np.ndarray]
            How long a patient can infect others (in days). If an array,
            it should be the infectious period for each population compartment.
        initial_R0: Union[int, float]
            Basic reproductive number of the disease
        hospitalization_probability: Union[float, np.ndarray]
            Probability that an infected person needs hospitalization.
            If an array, it should be the hospitalization probability
            for each population compartment.
        hospitalization_duration: Union[float, int]
            Average duration of a hospitalization in days.
        hospitalization_lag_from_onset: Union[float, int]
            Average time from the onset of symptoms to admission to hospital.
        icu_probability: Union[float, np.ndarray]
            Probability that an infected person needs hospitalization.
            If an array, it should be the probability for each population
            compartment.
        icu_duration: Union[float, int]
            Average duration of the need for intensive care in days.
        icu_lag_from_onset: Union[float, int]
            Average time from the onset of symptoms to admission to ICU.
        death_probability: Union[float, np.ndarray]
            Probability that an infected person dies from the disease.
            If an array, it should be the probability for each population
            compartment.
        death_lag_from_onset: Union[float, int]
            Average time from the onset of symptoms to death
        population: Union[int, float, np.ndarray]
            The total population. If an array, it should be the number of
            people in each population compartment.
        compartments: Optional[List[Any]]
            A description of each compartment of the population.
            For age-compartmentalized population this could be a list of
            the age for each compartment, e.g.,
            [ (0, 5), (5,10), (10,40), (40, 150), (150, 'inf') ]
        contacts_matrix: Optional[np.ndarray]
            A matrix C[i,j] describing the daily number of contacts a person of
            compartment 'i' has with the population of compartment 'j'.
        restrictions_function: Optional[Callable[[float],
                                        Union[float, np.ndarray]]]
            A function with signature `fun(time)` that outputs a matrix of
            the same shape as `contacts_matrix` or a float. At each
            timestep the `restrictions_function` is used to augment
            the infectivity rate matrix by a Hadamard product
            from the function's output.
        imported_cases_function: Optional[Callable]
        """
        # Set a single compartment if nothing was provided
        if compartments is not None:
            self.compartments = compartments
        else:
            self.compartments = ['All']

        self.num_compartments = len(self.compartments)

        # Save arguments inside the instance
        self.incubation_period = self._fix_size(incubation_period)
        self.infectious_period = self._fix_size(infectious_period)
        self.initial_R0 = initial_R0
        self.hospitalization_probability = self._fix_size(hospitalization_probability)
        self.hospitalization_duration = hospitalization_duration
        self.hospitalization_lag_from_onset = hospitalization_lag_from_onset
        self.icu_probability = self._fix_size(icu_probability)
        self.icu_duration = icu_duration
        self.icu_lag_from_onset = icu_lag_from_onset
        self.death_probability = self._fix_size(death_probability)
        self.death_lag_from_onset = death_lag_from_onset
        # Sanity checking on the population argument
        if isinstance(population, (int, float)):
            assert self.num_compartments == 1
        elif isinstance(population, list):
            population = np.array(population)
            assert population.size == self.num_compartments
        elif isinstance(population, np.ndarray):
            assert population.size == self.num_compartments
        self.population = self._fix_size(population)

        # Sanity checking on the contacts_matrix argument
        if contacts_matrix is not None:
            assert contacts_matrix.shape[0] == len(self.compartments)
            assert contacts_matrix.shape[1] == len(self.compartments)
        else:
            contacts_matrix = np.ones((self.num_compartments, self.num_compartments))

        self.infectivity_matrix = self._compute_infectivity_matrix(contacts_matrix)

        self.restrictions_function = restrictions_function
        self.imported_cases_function = imported_cases_function
        self._inf_matrix = self.infectivity_matrix.copy()  # Allocate once the infectivity matrix used in computation
        # Initial state
        self.Y0: Optional[np.ndarray] = None

        # Computed solution
        self.SEIR_solution: Optional[Callable[[np.ndarray], np.ndarray]] = None

    def update_parameters(
        self,
        *,
        incubation_period: Optional[Union[int, float, np.ndarray]] = None,
        infectious_period: Optional[Union[int, float, np.ndarray]] = None,
        initial_R0: Optional[Union[int, float]] = None,
        hospitalization_probability: Optional[Union[float, np.ndarray]] = None,
        hospitalization_duration: Optional[Union[float, np.ndarray]] = None,
        hospitalization_lag_from_onset: Optional[Union[float, np.ndarray]] = None,
        icu_probability: Optional[Union[float, np.ndarray]] = None,
        icu_duration: Optional[Union[float, np.ndarray]] = None,
        icu_lag_from_onset: Optional[Union[float, np.ndarray]] = None,
        death_probability: Optional[Union[float, np.ndarray]] = None,
        death_lag_from_onset: Optional[Union[float, np.ndarray]] = None,
        population: Optional[Union[float, np.ndarray]] = None,
        compartments: Optional[List[Any]] = None,
        contacts_matrix: Optional[np.ndarray] = None,
        restrictions_function: Optional[Callable[[float], Union[float, np.ndarray]]] = None,
        imported_cases_function: Optional[Callable] = None
    ):
        """
        Initializes the SEIR models parameters and computes the infectivity
        rate from the contacts matrix, R0, and infective_duration. Supports
        compartmentalizing the population to, e.g., age-groups.

        Keyword arguments
        -----------------
        incubation_period: Union[int, float, np.ndarray]
            Incubation period of the disease in days. If an array,
            it should be the incubation period for each population compartment.
        infectious_period: Union[int, float, np.ndarray]
            How long a patient can infect others (in days). If an array,
            it should be the infectious period for each population compartment.
        initial_R0: Union[int, float]
            Basic reproductive number of the disease
        hospitalization_probability: Union[float, np.ndarray]
            Probability that an infected person needs hospitalization.
            If an array, it should be the hospitalization probability
            for each population compartment.
        hospitalization_duration: Union[float, int]
            Average duration of a hospitalization in days.
        hospitalization_lag_from_onset: Union[float, int]
            Average time from the onset of symptoms to admission to hospital.
        icu_probability: Union[float, np.ndarray]
            Probability that an infected person needs hospitalization.
            If an array, it should be the probability for each population
            compartment.
        icu_duration: Union[float, int]
            Average duration of the need for intensive care in days.
        icu_lag_from_onset: Union[float, int]
            Average time from the onset of symptoms to admission to ICU.
        death_probability: Union[float, np.ndarray]
            Probability that an infected person dies from the disease.
            If an array, it should be the probability for each population
            compartment.
        death_lag_from_onset: Union[float, int]
            Average time from the onset of symptoms to death
        population: Union[int, float, np.ndarray]
            The total population. If an array, it should be the number of
            people in each population compartment.
        compartments: Optional[List[Any]]
            A description of each compartment of the population.
            For age-compartmentalized population this could be a list of
            the age for each compartment, e.g.,
            [ (0, 5), (5,10), (10,40), (40, 150), (150, 'inf') ]
        contacts_matrix: Optional[np.ndarray]
            A matrix C[i,j] describing the daily number of contacts a person of
            compartment 'i' has with the population of compartment 'j'.
        restrictions_function: Optional[Callable[[float],
                                        Union[float, np.ndarray]]]
            A function with signature `fun(time)` that outputs a matrix of
            the same shape as `contacts_matrix` or a float. At each
            timestep the `restrictions_function` is used to augment
            the infectivity rate matrix by a Hadamard product
            from the function's output.
        imported_cases_function: Optional[Callable]
        """
        raise NotImplementedError("TODO")

    def _compute_infectivity_matrix(self, contacts_matrix: np.ndarray) -> np.ndarray:
        r"""
        Computes and returns the infectivity matrix $\mathcal{I}_{a,a^*}$
        from the contacts matrix, R_0, and infective duration.

        See https://covid19.solanpaa.fi/#infectivity-rate for details.

        Input
        -----
        contacts_matrix: np.ndarray
            A matrix C[i,j] describing the daily number of contacts a person of
            compartment 'i' has with the population of compartment 'j'.

        Output
        ------
        infectivity_matrix: np.ndarray
        """
        normalization = 1 / self.infectious_period * \
                        self.initial_R0 * self.population.sum() / \
                        (self.population @ contacts_matrix).sum()
        # Symmetrize the contact matrix here since if
        # compartment 'i' has contacts with compartment 'j',
        # it should also be vice versa
        return normalization * 0.5 * (contacts_matrix + contacts_matrix.T)

    def _fix_size(self, x: Union[np.ndarray, float, int]) -> np.ndarray:
        """
        Fixes the size of the input to have
        the same size as there are compartments.

        Input
        -----
        x: Union[np.ndarray, float, int]
            The input parameter

        Output
        ------
        np.ndarray:
            the output array with same size as there are
            compartments in the model
        """
        if isinstance(x, list):
            x = np.array(x)

        if isinstance(x, (int, float)):
            return np.ones(self.num_compartments) * x
        else:
            assert x.size == self.num_compartments
            return x

    def __call__(self, t, Y):
        """
        Computes dY/dt of the SEIR model.

        Input
        -----
        t: float
            Time
        Y: np.ndarray
            The state of the system at time `t`, i.e.,
            Y = [S_0, S_1, ..., S_{N-1},
                 E_0, E_1, ..., E_{N-1},
                 I_0, I_1, ..., I_{N-1},
                 R_0, R_1, ..., R_{N-1}],
            where N is the number of compartments, and S_i, E_i, I_i, R_i
            are the number of susceptible, exposed, infected, and removed
            people at time `t` of the compartment `i`.
        Output
        ------
        dY/dt : np.ndarray with same shape as the input `Y`
        """
        if self.restrictions_function:
            self._inf_matrix[:, :] = np.multiply(self.restrictions_function(t), self.infectivity_matrix)

        # Create views to the state Y
        Sa = Y[:self.num_compartments]
        Ea = Y[self.num_compartments:2 * self.num_compartments]
        Ia = Y[2 * self.num_compartments:3 * self.num_compartments]

        res = np.zeros(4 * self.num_compartments)
        dS_dt = res[:self.num_compartments]
        dE_dt = res[self.num_compartments:2 * self.num_compartments]
        dI_dt = res[2 * self.num_compartments:3 * self.num_compartments]
        dR_dt = res[3 * self.num_compartments:]

        dS_dt[:] = -np.divide(Sa, self.population) * (self._inf_matrix @ Ia)
        dE_dt[:] = -dS_dt - np.divide(Ea, self.incubation_period)
        dI_dt[:] = np.divide(Ea, self.incubation_period) - np.divide(Ia, self.infectious_period)
        dR_dt[:] = np.divide(Ia, self.infectious_period)

        if self.imported_cases_function:
            DS, DE, DI = self.imported_cases_function(t)
            dS_dt += DS
            dE_dt += DE
            dI_dt += DI
        return res

    def set_initial_state(
        self,
        population_exposed: Union[int, float, np.ndarray],
        population_infected: Union[int, float, np.ndarray],
        probabilities: bool = False
    ):
        """
        Sets the initial state of the population system.

        Input
        -----
        population_exposed:
            If `probabilities` is True, this is the probability
            (or probabilities for each compartment) that a person
            is initially in the Exposed state.

            If `probabilities` is False, this is the
            number (or numbers for each compartment) of
            persons is initially in the Exposed state.
        population_infected:
            If `probabilities` is True, this is the probability
            (or probabilities for each compartment) that a person
            is initially in the Infected state.

            If `probabilities` is False, this is the
            number (or numbers for each compartment) of
            persons is initially in the Infected state.
        probabilities: bool, optional
            If True, the previous arguments are interpreted as
            probabilities. If False, they are interpreted as
            the number of people.
        """
        if probabilities:
            exposed = np.multiply(population_exposed, self.population)
            infected = np.multiply(population_infected, self.population)
        else:

            if isinstance(population_exposed, (int, float)):
                exposed = self._fix_size(population_exposed) / \
                          self.num_compartments
            elif isinstance(population_exposed, (list, np.ndarray)):
                population_exposed = np.array(population_exposed)
                assert population_exposed.size == self.num_compartments
                exposed = population_exposed

            if isinstance(population_infected, (int, float)):
                infected = self._fix_size(population_infected) / self.num_compartments
            elif isinstance(population_infected, (list, np.ndarray)):
                population_infected = np.array(population_infected)
                assert population_infected.size == self.num_compartments
                infected = population_infected

        susceptible = self.population - exposed - infected

        removed = np.zeros(self.num_compartments)

        self.Y0 = np.concatenate([susceptible, exposed, infected, removed])

    def simulate(
        self, max_simulation_time: Union[int, float], t0: float = 0.0, max_step: float = 0.5,
        method: Text = 'DOP853'
    ):
        """
        Simulates the SEIR model.

        Input
        -----
        days_to_simulate: Union[int, float]
            How many days forward to simulate the model
        """
        assert self.Y0 is not None
        solution = solve_ivp(
            fun=self,
            t_span=[t0, max_simulation_time],
            y0=self.Y0,
            dense_output=True,
            max_step=max_step,
            method=method,
        )

        # Create a callable returning the solution of the model
        def SEIR_solution(time: np.ndarray):
            postime_mask = time >= t0
            if np.sum(postime_mask) == time.size:
                return np.swapaxes(solution.sol(time), 0, 1)
            elif np.sum(postime_mask) == 0:
                negsol = np.repeat(np.expand_dims(self.Y0, 0), time.size, 0)
                return negsol
            else:
                possol = np.swapaxes(solution.sol(time[postime_mask]), 0, 1)
                negsol = np.repeat(np.expand_dims(self.Y0, 0), time.size - sum(postime_mask), 0)
                return np.vstack([negsol, possol])

        self.t0 = t0
        self.SEIR_solution = SEIR_solution

    def evaluate_solution(self, time: np.ndarray, only_real_results: bool = True):
        """
        Evaluates the results of a previously computed simulation.

        Input
        -----
        time: np.ndarray
            Time instances for which to evaluate the simulated model.
            Expected to be equidistant.
        only_real_results: bool
            If true, all results for which there is not enough history are set to NaN.

        Returns
        -------
        pd.DataFrame with columns

            - time
            - (susceptible, <Compartment 1>)
            - (susceptible, <Compartment 2>)
            - ...
            - susceptible
            - (exposed, <Compartment 1>)
            - (exposed, <Compartment 2>)
            - ...
            - exposed
            - (infected (active), <Compartment 1>)
            - (infected (active), <Compartment 2>)
            - ...
            - infected (active)
            - (infected (total), <Compartment 1>)
            - (infected (total), <Compartment 2>)
            - ...
            - infected (total)
            - (removed, <Compartment 1>)
            - (removed, <Compartment 2>)
            - ...
            - removed
            - (hospitalized (active), <Compartment 1>)
            - (hospitalized (active), <Compartment 2>)
            - hospitalized (active)
            - ...
            - (in ICU, <Compartment 1>)
            - (in ICU, <Compartment 2>)
            - ...
            - in ICU
            - (deaths, <Compartment 1>)
            - (deaths, <Compartment 2>)
            - ...
            - deaths

            Note that deaths is the number of deaths since t=t0
        """
        # Evaluate SEIR model results
        assert self.SEIR_solution is not None

        dt = time[1] - time[0]

        # Positive time mask
        idx_pos = np.where(time >= self.t0)[0]
        idx_neg = np.where(time < self.t0)[0]
        SEIR = self.SEIR_solution(time)
        S, E, I, R = np.split(SEIR, 4, axis=-1)
        S0, E0, I0, R0 = np.split(self.Y0, 4)

        # Compute the cumulative sum of infected people
        if idx_pos.size > 0:
            Ipos_zero = trapz(
                np.divide(np.vstack([E0, E[idx_pos[0], :]]), self.incubation_period),
                x=[0, time[idx_pos[0]]],
                axis=0
            )
            Idiff_pos = cumtrapz(
                np.divide(E[idx_pos, :], self.incubation_period), x=time[idx_pos], axis=0, initial=0
            ) + Ipos_zero
        if idx_neg.size > 0:
            Idiff_neg = I[idx_neg, :] - I0
        if idx_neg.size > 0 and idx_pos.size > 0:
            Icumulative = I0 + np.concatenate([Idiff_neg, Idiff_pos], axis=0)
        elif idx_pos.size > 0:
            Icumulative = I0 + Idiff_pos
        elif idx_neg.size > 0:
            Icumulative = I0 + Idiff_neg
        else:
            raise RuntimeError("Should not happen.")

        # Compute the number of hospitalized people for each day
        htime_min = time.min() - self.hospitalization_duration
        htime_max = time.max() + self.hospitalization_duration
        htime_beg = np.arange(htime_min, time[0] + dt / 2, dt)
        htime_end = np.arange(time[-1], htime_max + dt / 2, dt)
        htime = np.concatenate([htime_beg, time, htime_end])
        Shl, Ehl, Ihl, Rhl = np.split(self.SEIR_solution(htime - self.hospitalization_lag_from_onset), 4, axis=-1)
        H_new_cases_rate = \
            np.multiply(self.hospitalization_probability, np.divide(Ehl, self.incubation_period))
        Hwindow = np.ones(int(round(self.hospitalization_duration / dt)))
        H_active_cases = np.stack(
            [np.convolve(H_new_cases_rate[:, i], Hwindow, mode='same') * dt for i in range(self.num_compartments)],
            axis=-1
        )
        if only_real_results:
            H_t0 = self.t0 + self.hospitalization_lag_from_onset + self.hospitalization_duration
            H_active_cases[htime < H_t0] = np.nan
        H_active_cases = H_active_cases[htime_beg.size:htime_beg.size + time.size]

        # Compute the number of people in ICU for each day
        itime_min = time.min() - self.icu_duration
        itime_max = time.max() + self.icu_duration
        itime_beg = np.arange(itime_min, time[0] + dt / 2, dt)
        itime_end = np.arange(time[-1], itime_max + dt / 2, dt)
        itime = np.concatenate([itime_beg, time, itime_end])
        SEIR_icu_lag = self.SEIR_solution(itime - self.icu_lag_from_onset)
        E_icu_lag = np.split(SEIR_icu_lag, 4, axis=-1)[2]
        ICU_new_cases_rate = \
            np.multiply(self.icu_probability, np.divide(E_icu_lag, self.incubation_period))

        ICUwindow = np.ones(int(round(self.icu_duration / dt)))
        ICU_active_cases = np.stack(
            [np.convolve(ICU_new_cases_rate[:, i], ICUwindow, mode='same') * dt for i in range(
                self.num_compartments)],
            axis=-1
        )
        if only_real_results:
            ICU_t0 = self.t0 + self.icu_lag_from_onset + self.icu_duration
            ICU_active_cases[itime < ICU_t0] = np.nan
        ICU_active_cases = ICU_active_cases[itime_beg.size:itime_beg.size + time.size]

        # Compute the total number of deaths
        SEIR_death_lag = self.SEIR_solution(time - self.death_lag_from_onset)
        E_death_lag = np.split(SEIR_death_lag, 4, axis=-1)[1]
        deaths = cumtrapz(
            np.multiply(self.death_probability, np.divide(E_death_lag, self.incubation_period)),
            x=time,
            axis=0,
            initial=0
        )
        # Set deaths to NaN for t<=t0, as we have no real knowledge of those
        if only_real_results:
            deaths_t0 = self.t0 + self.death_lag_from_onset
            idx_closest_to_deaths_t0 = np.fabs(time - deaths_t0).argmin()
            deaths -= deaths[idx_closest_to_deaths_t0]
            deaths[time < deaths_t0] = np.nan

        # Compute the total cases over all compartments
        Sall = np.expand_dims(np.sum(S, axis=-1), -1)
        Eall = np.expand_dims(np.sum(E, axis=-1), -1)
        Iall = np.expand_dims(np.sum(I, axis=-1), -1)
        Icum_all = np.expand_dims(np.sum(Icumulative, axis=-1), -1)
        Rall = np.expand_dims(np.sum(R, axis=-1), -1)
        H_all = np.expand_dims(np.sum(H_active_cases, axis=-1), -1)
        ICU_all = np.expand_dims(np.sum(ICU_active_cases, axis=-1), -1)
        deaths_all = np.expand_dims(np.sum(deaths, axis=-1), -1)

        # Form the results dataframe
        data = np.concatenate(
            [
                np.expand_dims(time, -1), S, Sall, E, Eall, I, Iall, Icumulative, Icum_all, R, Rall, H_active_cases,
                H_all, ICU_active_cases, ICU_all, deaths, deaths_all
            ],
            axis=-1
        )
        columns = [
            'susceptible', 'exposed', 'infected (active)', 'infected (total)', 'removed', 'hospitalized (active)',
            'in ICU', 'deaths'
        ]
        all_columns = ['time'] + list(
            itertools.chain.from_iterable(
                [list(itertools.product([colname], self.compartments)) + [colname] for colname in columns]
            )
        )

        return pd.DataFrame(data, columns=all_columns)
