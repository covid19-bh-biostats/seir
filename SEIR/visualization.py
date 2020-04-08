from typing import Dict, List, Optional, Text

import matplotlib.pyplot as plt
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


def _plot_restrictions(ax, restrictions_info):
    for i, info in enumerate(restrictions_info):
        ax.plot([info['begins'], info['ends']], [i, i], linewidth=3)
    ax.set_ylim(-1, i + 1)
    ax.set_yticks(list(range(len(restrictions_info))))
    ax.set_yticklabels([ri['title'] for ri in restrictions_info])
    ax.tick_params(axis='y', which='major', labelsize=8)


def visualize_seir_computation(results: pd.DataFrame,
                               compartments: List[Text],
                               restrictions_info: Optional[List[Dict]] = None,
                               show_individual_compartments: bool = False):
    """Visualizes the SEIR computation"""

    if show_individual_compartments:
        w, h = plt.figaspect(2)
        if restrictions_info:
            w *= 2
        fig = plt.figure(figsize=(w, h))

        gs = GridSpec(1, 2, fig, width_ratios=[5, 1])
        gsp = GridSpecFromSubplotSpec(
            5,
            1,
            gs[0],
            hspace=0,
            height_ratios=[1 if restrictions_info else 0, 3, 3, 3, 3])

        if restrictions_info:
            ax = fig.add_subplot(gsp[0])
            _plot_restrictions(ax, restrictions_info)
        else:
            ax = None

        if ax:
            ax = fig.add_subplot(gsp[1], sharex=ax)
        else:
            ax = fig.add_subplot(gsp[1])

        _plot_compartment_subplot(ax, 'susceptible', results)

        ax = fig.add_subplot(gsp[2], sharex=ax)
        _plot_compartment_subplot(ax, 'exposed', results)

        ax = fig.add_subplot(gsp[3], sharex=ax)
        _plot_compartment_subplot(ax, 'infected (active)', results)

        ax = fig.add_subplot(gsp[4], sharex=ax)
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

    # Visualize restrictions (if any)
    if restrictions_info:
        ylim = ax.get_ylim()
        for rinfo in restrictions_info:

            ax.fill_between(
                [rinfo['begins'], rinfo['ends']],
                [0, 0],
                [ylim[1], ylim[1]],
                alpha=0.3,
                label=rinfo['title'])
    ax.legend()
    ax.set_xlabel('time (days)')
    ax.set_ylabel('# of people')

    ax.yaxis.set_major_formatter(EngFormatter())
    fig.tight_layout()
    plt.show()
