"""Test scoring of screw axes."""
from dials.array_family import flex
from dials.algorithms.symmetry.absences.screw_axes import (
    ScrewAxis41c,
    ScrewAxis42c,
    ScrewAxis61c,
    ScrewAxis21c,
)

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
    r["variance"] = flex.double(s) ** 2
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
    r["variance"] = flex.double(s) ** 2
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
    r["variance"] = flex.double(s) ** 2
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
