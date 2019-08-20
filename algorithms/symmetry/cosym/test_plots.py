from __future__ import absolute_import, division, print_function

from dials.algorithms.symmetry.cosym import plots
from scitbx.array_family import flex


def test_plot_coords():
    coords = flex.double([0, 1, 0, 1, 1, 0, 1, 0])
    coords.reshape(flex.grid((4, 2)))
    labels = flex.int([0, 0, 1, 1])
    d = plots.plot_coords(coords, labels=labels)
    assert "cosym_coordinates" in d
    assert set(d["cosym_coordinates"]) == {"layout", "data"}
    assert d["cosym_coordinates"]["data"][0]["x"] == [0.0, 0.0]
    assert d["cosym_coordinates"]["data"][0]["y"] == [1.0, 1.0]
    assert d["cosym_coordinates"]["data"][1]["x"] == [1.0, 1.0]
    assert d["cosym_coordinates"]["data"][1]["y"] == [0.0, 0.0]


def test_plot_rij_histogram():
    rij_matrix = flex.random_double(16)
    d = plots.plot_rij_histogram(rij_matrix)
    assert "cosym_rij_histogram" in d
    assert sum(d["cosym_rij_histogram"]["data"][0]["y"]) == 16
