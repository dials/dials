from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex

from dials.algorithms.symmetry.cosym import plots


def test_plot_coords():
    coords = flex.double([0, 1, 0, 1, 1, 0, 1, 0])
    coords.reshape(flex.grid((4, 2)))
    labels = flex.int([0, 0, 1, 1])
    d = plots.plot_coords(coords, labels=labels)
    assert "coordinates" in d
    assert d["coordinates"].keys() == ["layout", "data"]
    assert d["coordinates"]["data"][0]["x"] == [0.0, 0.0]
    assert d["coordinates"]["data"][0]["y"] == [1.0, 1.0]
    assert d["coordinates"]["data"][1]["x"] == [1.0, 1.0]
    assert d["coordinates"]["data"][1]["y"] == [0.0, 0.0]


def test_plot_rij_histogram():
    rij_matrix = flex.random_double(16)
    d = plots.plot_rij_histogram(rij_matrix)
    assert "rij_histogram" in d
    sum(d["rij_histogram"]["data"][0]["y"]) == 16
