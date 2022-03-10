from __future__ import annotations

import pytest

from dxtbx.model import BeamFactory, CrystalFactory, DetectorFactory, Experiment, Scan


@pytest.fixture
def test_experiment():
    d = {
        "__id__": "crystal",
        "real_space_a": [14.963210089244596, -22.599814679318, 51.02946725220764],
        "real_space_b": [-19.963976860932235, -51.503385430151205, -16.955728379753463],
        "real_space_c": [135.29560393219694, -34.371677531924206, -54.89475471853507],
        "space_group_hall_symbol": " P 4",
        "A_at_scan_points": [
            [
                0.004481726844090139,
                -0.005980612987053365,
                0.006013325470974739,
                -0.006768741824936281,
                -0.015428970379357122,
                -0.0015280122438480544,
                0.01528745348419002,
                -0.005078101688718203,
                -0.0024394384982453095,
            ]
        ],
    }

    crystal = CrystalFactory.from_dict(d)

    beam_d = {
        "direction": [-2.4882593300783137e-06, -0.0, 0.9999999999969044],
        "transmission": 1.0,
        "polarization_normal": [0.0, 1.0, 0.0],
        "divergence": 0.0,
        "polarization_fraction": 0.999,
        "flux": 0.0,
        "sigma_divergence": 0.0,
        "wavelength": 0.9762499999999994,
    }

    beam = BeamFactory.from_dict(beam_d)
    scan = Scan(image_range=[0, 1], oscillation=[0.0, 0.01])

    detector_dict = {
        "hierarchy": {
            "origin": [0.0, 0.0, 0.0],
            "fast_axis": [1.0, 0.0, 0.0],
            "name": "",
            "raw_image_offset": [0, 0],
            "slow_axis": [0.0, 1.0, 0.0],
            "material": "",
            "mask": [],
            "thickness": 0.0,
            "mu": 0.0,
            "gain": 1.0,
            "trusted_range": [0.0, 0.0],
            "image_size": [0, 0],
            "px_mm_strategy": {"type": "SimplePxMmStrategy"},
            "identifier": "",
            "type": "",
            "children": [{"panel": 0}],
            "pixel_size": [0.0, 0.0],
        },
        "panels": [
            {
                "origin": [-210.66631009735772, 205.7063614421482, -263.8386975038205],
                "fast_axis": [
                    0.9999973940105483,
                    -0.0016357501034268717,
                    -0.0015925745544149894,
                ],
                "name": "Panel",
                "raw_image_offset": [0, 0],
                "slow_axis": [
                    -0.0016426481736367285,
                    -0.999989234013669,
                    -0.004339765400707805,
                ],
                "material": "Si",
                "mask": [
                    [488, 1, 494, 2527],
                    [982, 1, 988, 2527],
                    [1476, 1, 1482, 2527],
                    [1970, 1, 1976, 2527],
                    [1, 196, 2463, 212],
                    [1, 408, 2463, 424],
                    [1, 620, 2463, 636],
                    [1, 832, 2463, 848],
                    [1, 1044, 2463, 1060],
                    [1, 1256, 2463, 1272],
                    [1, 1468, 2463, 1484],
                    [1, 1680, 2463, 1696],
                    [1, 1892, 2463, 1908],
                    [1, 2104, 2463, 2120],
                    [1, 2316, 2463, 2332],
                ],
                "thickness": 0.32,
                "mu": 3.9220322752480934,
                "gain": 1.0,
                "trusted_range": [-1.0, 161977.0],
                "image_size": [2463, 2527],
                "px_mm_strategy": {"type": "ParallaxCorrectedPxMmStrategy"},
                "identifier": "",
                "type": "SENSOR_PAD",
                "pixel_size": [0.17200000000000001, 0.17200000000000001],
            }
        ],
    }

    detector = DetectorFactory.from_dict(detector_dict)

    expt = Experiment(beam=beam, crystal=crystal, scan=scan, detector=detector)
    return expt
