from typing import Any, List

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.cm as cm
from matplotlib.ticker import EngFormatter
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
import pandas as pd
import numpy as np


def _plot_compartment_subplot(ax, observable, results):
    """
    Plots 'observable' of all compartments to a single axis.
    """
    compartments = np.unique(
        [v[1] for v in results.columns.values if isinstance(v, tuple)])
    colors_from_cmap = compartments.size > 6
    cmap = cm.get_cmap('viridis')
    lines = []
    for i, c in enumerate(compartments):
        if colors_from_cmap:
            lines += ax.plot(results['time'],
                             results[(observable, c)],
                             label=c,
                             color=cmap(i / compartments.size))
        else:
            lines += ax.plot(results['time'],
                             results[(observable, c)],
                             label=c)
    ax.set_ylabel(observable.capitalize())
    ax.yaxis.set_major_formatter(EngFormatter())
    return lines


def visualize_seir_computation(results: pd.DataFrame,
                               compartments: List[Any],
                               show_individual_compartments=False
                              ):
    """Visualizes the SEIR computation"""

    if show_individual_compartments:
        w, h = plt.figaspect(2)
        fig = plt.figure(figsize=(w, h))

        gs = GridSpec(1, 2, fig, width_ratios=[5, 1])
        gsp = GridSpecFromSubplotSpec(4, 1, gs[0], hspace=0)

        ax = fig.add_subplot(gsp[0])
        _plot_compartment_subplot(ax, 'susceptible', results)

        ax = fig.add_subplot(gsp[1], sharex=ax)
        _plot_compartment_subplot(ax, 'exposed', results)

        ax = fig.add_subplot(gsp[2], sharex=ax)
        _plot_compartment_subplot(ax, 'infected (active)', results)

        ax = fig.add_subplot(gsp[3], sharex=ax)
        lines = _plot_compartment_subplot(ax, 'deaths', results)

        ax.yaxis.set_major_formatter(EngFormatter())
        ax = fig.add_subplot(gs[1])
        ax.legend(lines, compartments)
        ax.set_xticks(())
        ax.set_yticks(())
        ax.set_axis_off()
        fig.tight_layout()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(results['time'], results['exposed'], label='exposed')
    ax.plot(results['time'],
            results['infected (active)'],
            label='infected (active)')

    ax.plot(results['time'],
            results['hospitalized (active)'],
            label='hospitalized (active)')

    ax.plot(results['time'], results['in ICU'], label='in ICU')

    ax.plot(results['time'], results['deaths'], label='deaths', color='k')
    ax.legend()
    ax.set_xlabel('time (days)')
    ax.set_ylabel('# of people')
    ax.yaxis.set_major_formatter(EngFormatter())
    fig.tight_layout()
    plt.show()
