from cctbx import miller
from dxtbx.model import ExperimentList
from dials.algorithms.spot_prediction import ScanStaticReflectionPredictor
from dials.util import missing_reflections


def test_connected_components(dials_data):
    experiment = ExperimentList.from_file(
        dials_data("centroid_test_data").join("experiments.json").strpath
    )[0]

    experiment.scan.set_image_range((1, 100))
    image_ranges = [(1, 9), (1, 100), (1, 1000)]
    expected_ms_sizes = [[755], [242, 14, 10, 5, 2, 2, 2], []]
    for image_range, expected_sizes in zip(image_ranges, expected_ms_sizes):
        experiment.scan.set_image_range(image_range)
        predict = ScanStaticReflectionPredictor(experiment, dmin=3, margin=1)
        refl = predict.for_ub(experiment.crystal.get_A())
        miller_set = miller.set(
            experiment.crystal.get_crystal_symmetry(),
            refl["miller_index"],
            anomalous_flag=False,
        )
        miller_array = miller_set.d_spacings().resolution_filter(d_min=3)
        complete_set, unique_ms = missing_reflections.connected_components(miller_array)
        assert len(unique_ms) == len(expected_sizes)
        assert [ms.size() for ms in unique_ms] == expected_sizes
