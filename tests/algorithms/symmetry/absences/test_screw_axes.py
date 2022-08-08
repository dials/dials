"""Test scoring of screw axes."""

from __future__ import annotations

from dials.algorithms.symmetry.absences.screw_axes import (
    ScrewAxis21c,
    ScrewAxis31c,
    ScrewAxis41c,
    ScrewAxis42c,
    ScrewAxis61c,
)
from dials.array_family import flex


# take some actual data for difficult cases - merged I, sigma values
def make_test_data_LCY_21c():
    """Real data from LCY (P212121)."""
    r = flex.reflection_table()
    miller_ax_vals = list(range(2, 17))
    i = [243.42, 0.95, 841.09, 4.08, 6780.07, 0.12, 19517.36, 0.10]
    i += [45.17, 1.60, 1056.85, 5.88, 7054.45, 7.77, 284.67]
    # separate to avoid black formatting
    s = [8.64, 0.80, 28.70, 1.87, 248.14, 1.30, 637.42, 1.57, 4.49, 2.05]
    s += [50.27, 3.37, 520.39, 4.84, 24.14]
    r["miller_index"] = flex.miller_index([(0, 0, j) for j in miller_ax_vals])
    r["variance"] = flex.pow2(flex.double(s))
    r["intensity"] = flex.double(i)
    return r


def make_test_data_thermo_61c():
    """Real data from thermolysin (P6122), has a large intensity in (0, 0, 11),
    (an outlier?) so good for testing resiliance to outliers."""
    r = flex.reflection_table()
    miller_ax_vals = list(range(2, 37))
    i = [-0.007, -0.011, 0.044, 0.065, 1.388, 0.009, -0.005, 0.039, 0.056]
    i += [1.619, 163.471, 0.024, -0.024, 0.047, 0.027, -0.023, 0.09, 0.037]
    i += [0.024, -0.007, -0.024, 0.323, 15.846, 0.177, 0.042, 0.139, 0.37]
    i += [1.185, 218.702, 0.041, 0.185, 0.09, 0.383, 0.476, 107.947]
    s = [0.005, 0.008, 0.013, 0.018, 0.046, 0.028, 0.033, 0.037, 0.043]
    s += [0.075, 0.608, 0.055, 0.059, 0.066, 0.101, 0.105, 0.108, 0.113, 0.125]
    s += [0.188, 0.192, 0.144, 0.38, 0.149, 0.146, 0.158, 0.17, 0.185, 1.288]
    s += [0.127, 0.088, 0.086, 0.059, 0.092, 0.975]
    r["miller_index"] = flex.miller_index([(0, 0, j) for j in miller_ax_vals])
    r["variance"] = flex.pow2(flex.double(s))
    r["intensity"] = flex.double(i)
    return r


def make_test_data_thaumatin_41c():
    """Real thaumatin data (P41212). Nice example of 41 screw axis,
    in trouble if we can't get this one right!"""
    r = flex.reflection_table()
    miller_ax_vals = list(range(1, 87))
    i = [-0.006, -0.027, -0.016, 0.094, 0.012, 0.041, 0.039, 605.708, -0.01]
    i += [0.058, 0.005, 406.319, 0.047, 0.043, 0.082, 0.754, 0.101, 0.126]
    i += [-0.17, 1.381, 0.149, -0.177, 0.175, 25.368, 0.007, 0.442, -0.753]
    i += [2944.402, 0.02, -0.41, 0.399, 1334.451, 0.35, 0.028, -0.353, 24.594]
    i += [0.459, 0.093, -0.666, 1068.914, -1.048, 0.882, 0.26, 391.129, 0.771]
    i += [-0.605, 1.923, 79.672, 0.539, 0.342, 0.673, 570.054, -0.624, -0.388]
    i += [-0.572, 1175.132, -0.764, 0.006, 1.027, 459.19, -0.116, 0.098, -0.186]
    i += [319.29, 0.591, 0.874, 0.265, 4.021, 0.246, 0.262, 0.552, 14.901]
    i += [-1.647, -1.776, 0.01, 57.969, 1.067, -0.751, 0.438, 0.98, 0.277]
    i += [0.317, -0.888, 11.128, 0.332, -0.608]
    s = [0.005, 0.016, 0.023, 0.055, 0.072, 0.088, 0.112, 12.797, 0.125]
    s += [0.148, 0.17, 8.9, 0.199, 0.274, 0.297, 0.385, 0.378, 0.345, 0.287]
    s += [0.371, 0.291, 0.343, 0.334, 1.276, 0.358, 0.574, 0.724, 61.67, 0.46]
    s += [0.475, 0.652, 40.835, 0.65, 0.674, 0.796, 1.533, 0.587, 0.663, 0.692]
    s += [23.811, 0.794, 0.707, 0.765, 9.989, 0.833, 0.818, 1.21, 3.354, 0.808]
    s += [0.836, 0.836, 14.031, 0.894, 0.914, 0.992, 26.717, 0.937, 0.88, 0.904]
    s += [12.025, 0.879, 0.886, 0.981, 9.202, 1.012, 1.478, 1.413, 1.596, 1.453]
    s += [1.411, 0.861, 2.561, 1.481, 1.162, 0.918, 3.367, 0.93, 1.03, 0.875]
    s += [1.0, 0.925, 0.918, 0.983, 1.631, 0.821, 1.35]
    r["miller_index"] = flex.miller_index([(0, 0, j) for j in miller_ax_vals])
    r["variance"] = flex.pow2(flex.double(s))
    r["intensity"] = flex.double(i)
    return r


def make_test_data_ptp1b_31c():
    """Real PTP1B data collected by the Hekstra Lab (P 31 2 1)."""
    r = flex.reflection_table()
    miller_ax_vals = [
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
        26,
        27,
        28,
        29,
        30,
        31,
        32,
        33,
        34,
        35,
        36,
        37,
        38,
        39,
        40,
        41,
        42,
        43,
        44,
        45,
        46,
        48,
        49,
    ]
    i = [
        -0.002,
        -0.018,
        1765.107,
        0.049,
        0.206,
        538.946,
        0.065,
        0.068,
        511.518,
        0.008,
        0.295,
        46.878,
        -0.106,
        0.374,
        118.479,
        0.509,
        0.554,
        84.699,
        0.084,
        0.335,
        413.968,
        -0.371,
        0.144,
        65.366,
        -0.185,
        30.277,
        0.283,
        1.785,
        9.957,
        1.041,
        0.3,
        2.142,
        0.199,
        0.118,
        8.955,
        0.269,
        0.744,
        2.981,
        -0.575,
        0.26,
        -0.41,
        0.032,
        -0.294,
        71.932,
        0.028,
        17.147,
        1.246,
    ]
    s = [
        0.007,
        0.05,
        3.662,
        0.125,
        0.14,
        2.866,
        0.232,
        0.242,
        3.421,
        0.321,
        0.381,
        1.236,
        0.397,
        0.465,
        2.174,
        0.478,
        0.595,
        2.072,
        0.624,
        0.693,
        4.767,
        0.807,
        0.787,
        4.427,
        0.955,
        1.806,
        1.125,
        1.13,
        1.367,
        1.259,
        1.287,
        1.203,
        1.153,
        1.341,
        1.506,
        1.128,
        1.344,
        1.381,
        1.16,
        1.26,
        1.363,
        1.194,
        1.224,
        3.229,
        1.286,
        3.01,
        1.487,
    ]

    r["miller_index"] = flex.miller_index([(0, 0, j) for j in miller_ax_vals])
    r["variance"] = flex.pow2(flex.double(s))
    r["intensity"] = flex.double(i)
    return r


def test_screw_axes_example_data():
    """Test some example data where we know the answer"""
    refls = make_test_data_LCY_21c()
    score = ScrewAxis21c().score_axis(refls)
    assert score > 0.95

    refls = make_test_data_thermo_61c()
    score = ScrewAxis61c().score_axis(refls)
    assert score > 0.95

    refls = make_test_data_thaumatin_41c()
    score_41 = ScrewAxis41c().score_axis(refls)
    score_42 = ScrewAxis42c().score_axis(refls)
    assert score_41 > 0.99
    assert score_42 > 0.99  # both should score highly

    refls = make_test_data_ptp1b_31c()
    score = ScrewAxis31c().score_axis(refls)
    assert score > 0.95


def test_not_screw_axes():
    """Test cases where a screw axis is not present."""
    reflection_table = flex.reflection_table()
    reflection_table["miller_index"] = flex.miller_index(
        [(0, 0, j) for j in range(1, 9)]
    )
    reflection_table["variance"] = flex.double(8, 1.0)
    reflection_table["intensity"] = flex.double(8, 5.0)
    score_41 = ScrewAxis41c().score_axis(reflection_table)
    score_42 = ScrewAxis42c().score_axis(reflection_table)

    assert score_41 < 0.01
    assert score_42 < 0.01

    # Example of weak pseudosymmetry
    reflection_table["intensity"] = flex.double(
        [3.0, 100.0, 3.0, 100.0, 3.0, 100.0, 3.0, 100.0]
    )
    score_41 = ScrewAxis41c().score_axis(reflection_table)
    score_42 = ScrewAxis42c().score_axis(reflection_table)

    assert score_41 < 0.01
    assert score_42 < 0.01

    # All data too weak
    reflection_table["intensity"] = flex.double(8, 0.5)
    score_41 = ScrewAxis41c().score_axis(reflection_table)
    score_42 = ScrewAxis42c().score_axis(reflection_table)

    assert score_41 < 0.01
    assert score_42 < 0.01

    # No data for absent reflections
    reflection_table["miller_index"] = flex.miller_index(
        [(0, 0, i) for i in range(1, 17, 2)]
    )
    score_41 = ScrewAxis41c().score_axis(reflection_table)
    score_42 = ScrewAxis42c().score_axis(reflection_table)

    assert score_41 < 0.01
    assert score_42 < 0.01

    # No data for present reflections
    reflection_table["miller_index"] = flex.miller_index(
        [(0, 0, i) for i in range(2, 18, 2)]
    )
    score_41 = ScrewAxis41c().score_axis(reflection_table)
    score_42 = ScrewAxis42c().score_axis(reflection_table)

    assert score_41 < 0.01
    assert score_42 < 0.01
