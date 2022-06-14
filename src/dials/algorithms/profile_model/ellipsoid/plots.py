from __future__ import annotations

import numpy as np

from scitbx import matrix
from scitbx.array_family import flex


def plot_distance_from_ewald_sphere(experiment, reflection_table, prefix):
    """Plot distance from Ewald sphere"""

    s0 = matrix.col(experiment.beam.get_s0())
    s2 = reflection_table["s2"]
    D = flex.double(s0.length() - matrix.col(s).length() for s in s2)
    Dmean = flex.sum(D) / len(D)
    Dvar = flex.sum(flex.double([(d - Dmean) ** 2 for d in D])) / len(D)
    hist, bin_edges = np.histogram(
        D,
        bins=max(5, min(int(0.2 * len(s2)), 20)),
    )
    bin_centers = bin_edges[:-1] + np.diff(bin_edges) / 2
    plot = {
        f"{prefix}_epsilon_distribution": {
            "data": [
                (
                    {
                        "x": bin_centers.tolist(),
                        "y": hist.tolist(),
                        "type": "bar",
                    }
                )
            ],
            "layout": {
                "title": f"{prefix} epsilon distribution. <br>Mean(epsilon) = {Dmean:.2e}, Variance(epsilon) = {Dvar:.2e}",
                "xaxis": {"title": "Distance from Ewald sphere (epsilon)"},
                "yaxis": {"title": "Frequency"},
                "bargap": 0,
            },
        }
    }
    return plot


def plot_partiality(reflection_table):
    """Plot the partiality"""

    hist, bin_edges = np.histogram(
        reflection_table["partiality"],
        bins=max(5, min(int(0.2 * reflection_table.size()), 20)),
    )
    bin_centers = bin_edges[:-1] + np.diff(bin_edges) / 2
    plots_dict = {
        "partiality_distribution": {
            "data": [
                (
                    {
                        "x": bin_centers.tolist(),
                        "y": hist.tolist(),
                        "type": "bar",
                    }
                )
            ],
            "layout": {
                "title": "Partiality distribution",
                "xaxis": {"title": "Partiality"},
                "yaxis": {"title": "Frequency"},
                "bargap": 0,
            },
        }
    }
    return plots_dict
