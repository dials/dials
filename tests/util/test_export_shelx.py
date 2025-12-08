from __future__ import annotations

from os.path import join

import numpy as np

from dxtbx import flumpy
from scitbx.array_family import flex

from dials.util.export_shelx import _export_as_shelx_hklf2


def test_export_shelx_hklf2(tmp_path):
    intensities = np.array(
        [
            289.49048901,
            463.33289445,
            840.63032177,
            2268.21286727,
            470.66661385,
            2552.22515851,
            786.76105206,
            4221.55638505,
            849.60880376,
            1993.31277274,
        ]
    )

    sigmas = np.array(
        [
            21.16940444,
            24.40840987,
            31.1660734,
            50.19060457,
            23.89432053,
            52.96094679,
            30.89304753,
            67.10334197,
            30.90376238,
            46.55504274,
        ]
    )

    miller_indices = np.array(
        [
            [5, 1, 11],
            [5, 1, 9],
            [2, 2, 8],
            [-4, 0, 10],
            [-3, 3, 9],
            [-2, 4, 10],
            [-2, 0, 12],
            [-2, 4, 8],
            [-3, 1, 9],
            [-2, 2, 10],
        ],
        dtype=np.int32,
    )

    wavelengths = np.array(
        [
            0.90661541,
            1.06846022,
            1.18829931,
            0.60825594,
            0.64290777,
            0.65004844,
            0.68644674,
            0.70307848,
            0.7249276,
            0.74376255,
        ]
    )

    # Convert to flex arrays
    intensities = flumpy.from_numpy(intensities)
    sigmas = flumpy.from_numpy(sigmas)
    miller_indices = flumpy.miller_index_from_numpy(miller_indices)
    wavelengths = flumpy.from_numpy(wavelengths)
    batch_numbers = flex.int(len(wavelengths), 0)

    filename = join(tmp_path, "test.hkl")

    # Call the function
    _export_as_shelx_hklf2(
        filename=filename,
        miller_indices=miller_indices,
        intensities=intensities,
        sigmas=sigmas,
        batch_numbers=batch_numbers,
        wavelengths=wavelengths,
    )

    with open(filename, "r") as f:
        lines = f.readlines()

    assert len(lines) == len(intensities) + 1

    for i, line in enumerate(lines[:-1]):
        assert len(line.strip()) > 0
        h = int(line[0:4])
        k = int(line[4:8])
        l = int(line[8:12])
        intensity = float(line[12:20])
        sigma = float(line[20:28])
        batch = int(line[28:32])
        wavelength = float(line[32:40])

        np.testing.assert_almost_equal(intensity, round(intensities[i], 2), decimal=2)
        np.testing.assert_almost_equal(sigma, round(sigmas[i], 2), decimal=2)
        np.testing.assert_almost_equal(wavelength, round(wavelengths[i], 4), decimal=4)
        assert (h, k, l) == tuple(miller_indices[i])
        assert batch == 0

    # Validate last line is the termination line
    last = lines[-1]
    assert last == "   0   0   0    0.00    0.00   0  0.0000\n", (
        "Terminating line incorrect"
    )
