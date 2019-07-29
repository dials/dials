"""
Tests for intensity combination.
"""
from __future__ import absolute_import, division, print_function
import pytest
from mock import Mock
from dxtbx.model import Experiment, Crystal
from dials.array_family import flex
from dials.algorithms.scaling.scaling_utilities import calculate_prescaling_correction
from dials.algorithms.scaling.combine_intensities import (
    SingleDatasetIntensityCombiner,
    MultiDatasetIntensityCombiner,
)


@pytest.fixture(scope="module")
def test_exp_P1():
    """Create a mock experiments object."""
    exp = Experiment()
    exp_dict = {
        "__id__": "crystal",
        "real_space_a": [1.0, 0.0, 0.0],
        "real_space_b": [0.0, 1.0, 0.0],
        "real_space_c": [0.0, 0.0, 1.0],
        "space_group_hall_symbol": " P 1",
    }
    crystal = Crystal.from_dict(exp_dict)
    exp.crystal = crystal
    return exp


def generate_simple_table(prf=True):
    """Generate a reflection table for testing intensity combination.
    The numbers are contrived to make sum intensities agree well at high
    intensity but terribly at low and vice versa for profile intensities."""
    reflections = flex.reflection_table()
    reflections["miller_index"] = flex.miller_index(
        [
            (0, 0, 1),
            (0, 0, 1),
            (0, 0, 1),
            (0, 0, 1),
            (0, 0, 1),
            (0, 0, 2),
            (0, 0, 2),
            (0, 0, 2),
            (0, 0, 2),
            (0, 0, 2),
            (0, 0, 3),
            (0, 0, 3),
            (0, 0, 3),
            (0, 0, 3),
            (0, 0, 3),
            (0, 0, 4),
            (0, 0, 4),
            (0, 0, 4),
            (0, 0, 4),
            (0, 0, 4),
            (0, 0, 5),
            (0, 0, 5),
            (0, 0, 5),
            (0, 0, 5),
            (0, 0, 5),
        ]
    )
    reflections["inverse_scale_factor"] = flex.double(25, 1.0)
    # Contrive an example that should give the best cc12 when combined.
    # make sum intensities agree well at high intensity but terribly at low
    # and vice versa for profile intensities.
    # profile less consistent at high intensity here

    # sumless consistent at low intensity here
    reflections["intensity.sum.value"] = flex.double(
        [
            10000.0,
            11000.0,
            9000.0,
            8000.0,
            12000.0,
            500.0,
            5600.0,
            5500.0,
            2000.0,
            6000.0,
            100.0,
            50.0,
            150.0,
            75.0,
            125.0,
            30.0,
            10.0,
            2.0,
            35.0,
            79.0,
            1.0,
            10.0,
            20.0,
            10.0,
            5.0,
        ]
    )
    reflections["intensity.sum.variance"] = flex.double(
        [10000] * 5 + [5000] * 5 + [100] * 5 + [30] * 5 + [10] * 5
    )
    reflections.set_flags(flex.bool(25, False), reflections.flags.outlier_in_scaling)
    reflections.set_flags(flex.bool(25, True), reflections.flags.integrated)
    reflections["lp"] = flex.double(25, 0.5)
    if prf:
        reflections["intensity.prf.value"] = flex.double(
            [
                10000.0,
                16000.0,
                12000.0,
                6000.0,
                9000.0,
                5000.0,
                2000.0,
                1500.0,
                1300.0,
                9000.0,
                100.0,
                80.0,
                120.0,
                90.0,
                100.0,
                30.0,
                40.0,
                50.0,
                30.0,
                30.0,
                10.0,
                12.0,
                9.0,
                8.0,
                10.0,
            ]
        )
        reflections["intensity.prf.variance"] = flex.double(
            [1000] * 5 + [500] * 5 + [10] * 5 + [3] * 5 + [1] * 5
        )
    reflections = calculate_prescaling_correction(reflections)
    return reflections


def test_combine_intensities(test_exp_P1):
    """Test the combine intensities function for a single dataset"""
    reflections = generate_simple_table()
    scaler = Mock()
    scaler.reflection_table = reflections
    scaler.suitable_refl_for_scaling_sel = flex.bool(reflections.size(), True)
    scaler.outliers = flex.bool(reflections.size(), False)
    scaler.experiment = test_exp_P1
    scaler.space_group = test_exp_P1.crystal.get_space_group()
    scaler.params.reflection_selection.combine.Imid = None

    combiner = SingleDatasetIntensityCombiner(scaler)
    Imid = combiner.max_key
    intensity, variance = combiner.calculate_suitable_combined_intensities()

    # Imid being 1200.0 should be best for this contrived example
    assert Imid == 1200.0

    # Due to nature of crossover, just require 2% tolerance for this example
    assert list(intensity[0:5]) == pytest.approx(
        list(
            reflections["intensity.sum.value"][0:5]
            * reflections["prescaling_correction"][0:5]
        ),
        rel=2e-2,
    )
    assert list(intensity[20:25]) == pytest.approx(
        list(
            reflections["intensity.prf.value"][20:25]
            * reflections["prescaling_correction"][20:25]
        ),
        rel=2e-2,
    )

    assert list(variance[0:5]) == pytest.approx(
        list(
            reflections["intensity.sum.variance"][0:5]
            * reflections["prescaling_correction"][0:5] ** 2
        ),
        rel=2e-2,
    )
    assert list(variance[20:25]) == pytest.approx(
        list(
            reflections["intensity.prf.variance"][20:25]
            * reflections["prescaling_correction"][20:25] ** 2
        ),
        rel=2e-2,
    )

    combiner.max_key = 0  # prf
    intensity, variance = combiner.calculate_suitable_combined_intensities()
    assert list(intensity) == pytest.approx(
        list(reflections["intensity.prf.value"] * reflections["prescaling_correction"]),
        rel=2e-2,
    )
    assert list(variance) == pytest.approx(
        list(
            reflections["intensity.prf.variance"]
            * reflections["prescaling_correction"] ** 2
        ),
        rel=2e-2,
    )

    combiner.max_key = 1  # sum
    intensity, variance = combiner.calculate_suitable_combined_intensities()
    assert list(intensity) == pytest.approx(
        list(reflections["intensity.sum.value"] * reflections["prescaling_correction"]),
        rel=2e-2,
    )
    assert list(variance) == pytest.approx(
        list(
            reflections["intensity.sum.variance"]
            * reflections["prescaling_correction"] ** 2
        ),
        rel=2e-2,
    )


def test_combine_intensities_multi_dataset(test_exp_P1):
    """Test the combine intensities function for multiple datasets"""
    r1 = generate_simple_table()
    r1["partiality"] = flex.double(25, 1.0)
    r2 = generate_simple_table(prf=False)
    scaler1 = Mock()
    scaler1.reflection_table = r1
    scaler1.suitable_refl_for_scaling_sel = flex.bool(r1.size(), True)
    scaler1.outliers = flex.bool(r1.size(), False)
    scaler1.experiment = test_exp_P1
    scaler1.space_group = test_exp_P1.crystal.get_space_group()
    scaler1.params.reflection_selection.combine.Imid = None
    scaler2 = Mock()
    scaler2.reflection_table = r2
    scaler2.space_group = test_exp_P1.crystal.get_space_group()
    scaler2.suitable_refl_for_scaling_sel = flex.bool(r2.size(), True)
    scaler2.outliers = flex.bool(r2.size(), False)
    scaler2.experiment = test_exp_P1
    scaler2.params.reflection_selection.combine.Imid = None

    multiscaler = Mock()
    multiscaler.active_scalers = [scaler1, scaler2]
    multiscaler.experiment = test_exp_P1
    # multiscaler.space_group = test_exp_P1.crystal.get_space_group()
    multiscaler.params.reflection_selection.combine.Imid = None

    combiner = MultiDatasetIntensityCombiner(multiscaler)
    Imid = combiner.max_key

    # Imid = optimise_intensity_combination([r1, r2], test_exp_P1)
    assert pytest.approx(Imid) == 1200.0
