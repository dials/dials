# coding: utf-8
from __future__ import absolute_import, division, print_function

import numpy as np

from scitbx.array_family import flex


def plot_coords(coords, labels=None, key="cosym_coordinates"):

    coord_x = coords[:, 0:1].as_1d()
    coord_y = coords[:, 1:2].as_1d()
    assert coord_x.size() == coord_y.size(), (coord_x.size(), coord_y.size())

    if labels is None:
        labels = flex.int(len(coord_x), -1)

    unique_labels = set(labels)
    unique_labels = sorted(unique_labels)
    n_clusters = max(len(unique_labels) - (1 if -1 in unique_labels else 0), 1)

    # XXX should avoid relying on matplotlib here to determine colours
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    colours = plt.cm.Spectral(np.linspace(0, 1, n_clusters)).tolist()

    if -1 in unique_labels:
        colours.insert(0, (0, 0, 0, 1))
    data = []
    for k, col in zip(unique_labels, colours):
        isel = (labels == k).iselection()
        coord_x_sel = coord_x.select(isel)
        coord_y_sel = coord_y.select(isel)
        data.append(
            {
                "x": list(coord_x_sel),
                "y": list(coord_y_sel),
                "mode": "markers",
                "type": "scatter",
                "marker": {
                    "size": 2,
                    "alpha": 0.5,
                    "color": "rgb(%f,%f,%f)" % tuple(col[:3]),
                },
                "name": "Cluster %i" % k if k >= 0 else "Noise",
            }
        )
    d = {
        key: {
            "data": data,
            "layout": {
                "title": "Cosym coordinates",
                "xaxis": {"range": [-1, 1], "constrain": "domain"},
                "yaxis": {"range": [-1, 1], "scaleanchor": "x", "constrain": "domain"},
            },
            "help": """\
The outcome of the dials.cosym multi-dimensional scaling procedure projected on to two
dimensions. Each point corresponds to an individual data set, or symmetry copy of a data
set. The lengths of the vectors are inversely related to the amount of random error in
each data set, and can be interpreted as an estimate of the CC* values. The angular
separation between any pair, or groups, of vectors is a measure of the systematic
differences between the data sets, for example as a result of an indexing ambiguity,
or the presence of non-isomorphism.
""",
        }
    }

    return d


def plot_rij_histogram(rij_matrix, key="cosym_rij_histogram"):
    """Plot a histogram of the rij values.

    Args:
      plot_name (str): The file name to save the plot to.
        If this is not defined then the plot is displayed in interactive mode.
    """
    rij = rij_matrix.as_1d()
    rij = rij.select(rij != 0)
    hist = flex.histogram(
        rij,
        data_min=min(-1, flex.min(rij)),
        data_max=max(1, flex.max(rij)),
        n_slots=100,
    )

    d = {
        key: {
            "data": [
                {
                    "x": list(hist.slot_centers()),
                    "y": list(hist.slots()),
                    "type": "bar",
                    "name": "Rij histogram",
                }
            ],
            "layout": {
                "title": "Distribution of values in the Rij matrix",
                "xaxis": {"title": "r<sub>ij</sub>"},
                "yaxis": {"title": "Frequency"},
                "bargap": 0,
            },
            "help": """\
A histogram of the values of the Rij matrix of pairwise correlation coefficients. A
unimodal distribution of values may suggest that no indexing ambiguity is evident,
whereas a bimodal distribution can be indicative of the presence of an indexing
ambiguity.
""",
        }
    }
    return d
