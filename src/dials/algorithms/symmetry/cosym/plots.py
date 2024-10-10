from __future__ import annotations

import numpy as np


def plot_coords(coords, labels=None, key="cosym_coordinates", dim1=0, dim2=1):
    coord_x = coords[:, dim1]
    coord_y = coords[:, dim2]
    assert coord_x.size == coord_y.size, (coord_x.size, coord_y.size)

    if labels is None:
        labels = np.full(coord_x.shape[0], -1, dtype=int)

    unique_labels = set(labels)
    unique_labels = sorted(unique_labels)
    n_clusters = max(len(unique_labels) - (1 if -1 in unique_labels else 0), 1)

    # XXX should avoid relying on matplotlib here to determine colours
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    # Decimal places are used to avoid the start and end of the spectral map which are black (reserved for outliers) and white (does not show)

    colours = plt.cm.nipy_spectral(np.linspace(0.1, 0.9, n_clusters)).tolist()

    if -1 in unique_labels:
        colours.insert(0, (0, 0, 0, 1))
    data = []
    for k, col in zip(unique_labels, colours):
        isel = np.where(labels == k)[0]
        data.append(
            {
                "x": coord_x[isel].tolist(),
                "y": coord_y[isel].tolist(),
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
                "xaxis": {
                    "title": f"Axis {dim1}",
                    "range": [-1, 1],
                    "constrain": "domain",
                },
                "yaxis": {
                    "title": f"Axis {dim2}",
                    "range": [-1, 1],
                    "scaleanchor": "x",
                    "constrain": "domain",
                },
            },
            "help": """\
The outcome of the dials.cosym multi-dimensional scaling procedure projected on to two
dimensions. Each point corresponds to an individual data set, or symmetry copy of a data
set. The lengths of the vectors are inversely related to the amount of random error in
each data set, and can be interpreted as an estimate of the CC* values. The angular
separation between any pair, or groups, of vectors is a measure of the systematic
differences between the data sets, for example as a result of an indexing ambiguity,
or the presence of non-isomorphism. Where higher dimensions have been used for the analysis,
additional representations are given.
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
    hist, bin_edges = np.histogram(
        rij_matrix[rij_matrix != 0],
        bins=100,
        range=(min(-1, rij_matrix.min()), max(1, rij_matrix.max())),
    )
    bin_centers = bin_edges[:-1] + np.diff(bin_edges) / 2

    d = {
        key: {
            "data": [
                {
                    "x": bin_centers.tolist(),
                    "y": hist.tolist(),
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
