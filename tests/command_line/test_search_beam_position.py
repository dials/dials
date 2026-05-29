from __future__ import annotations

import os
import shutil
import subprocess
from collections import namedtuple

import mrcfile
import numpy as np
import pytest

import scitbx
from cctbx import uctbx
from dxtbx.serialize import load

from dials.command_line import search_beam_position

from ..algorithms.indexing.test_index import run_indexing

Gaussian = namedtuple("Gaussian", ["x0", "y0", "sigma_x", "sigma_y", "height"])


def test_search_i04_weak_data_image_range(mocker, run_in_tmp_path, dials_data):
    """Perform a beam-centre search and check that the output is sane."""

    data_dir = dials_data("i04_weak_data")
    reflection_file = data_dir / "full.pickle"
    experiments_file = data_dir / "experiments_import.json"

    args = [
        str(experiments_file),
        str(reflection_file),
        "image_range=1,10",
        "image_range=251,260",
        "image_range=531,540",
        "n_macro_cycles=4",
    ]
    from rstbx.indexing_api import dps_extended

    mocker.spy(dps_extended, "get_new_detector")
    search_beam_position.run(args)
    # Check that the last call to get_new_detector was with an
    # offset of close to zero.
    # The final call was to apply the "best" shift to the detector model before
    # returning the updated experiments.
    assert dps_extended.get_new_detector.call_args[0][1].elems == pytest.approx(
        (0, 0, 0), abs=3e-2
    )
    assert os.path.exists("optimised.expt")

    # Compare the shifts between the start and final detector models
    experiments = load.experiment_list(experiments_file, check_format=False)
    optimised_experiments = load.experiment_list("optimised.expt", check_format=False)
    detector_1 = experiments[0].detector
    detector_2 = optimised_experiments[0].detector
    shift = scitbx.matrix.col(detector_1[0].get_origin()) - scitbx.matrix.col(
        detector_2[0].get_origin()
    )
    assert shift.elems == pytest.approx((0.27, -0.12, 0.0), abs=1e-1)


def test_search_multiple(run_in_tmp_path, dials_data):
    """Perform a beam-centre search and check that the output is sane.

    Do the following:
    1. Run dials.search_beam_centre on two experiments and two pickled
    reflection tables, as output by dials.find_spots;
      a) Check that the program exits correctly;
      b) Check that it produces the expected output experiment.
    2. Check that the beam centre search has resulted in the expected shift
    in detector origin.
    """

    data_dir = dials_data("semisynthetic_multilattice")
    refl_path1 = str(data_dir / "ag_strong_1_50.refl")
    refl_path2 = str(data_dir / "bh_strong_1_50.refl")
    experiments_path1 = str(data_dir / "ag_imported_1_50.expt")
    experiments_path2 = str(data_dir / "bh_imported_1_50.expt")

    args = [experiments_path1, experiments_path2, refl_path1, refl_path2]
    search_beam_position.run(args)
    assert os.path.exists("optimised.expt")

    experiments = load.experiment_list(experiments_path1, check_format=False)
    optimised_experiments = load.experiment_list("optimised.expt", check_format=False)
    detector_1 = experiments[0].detector
    detector_2 = optimised_experiments[0].detector
    shift = scitbx.matrix.col(detector_1[0].get_origin()) - scitbx.matrix.col(
        detector_2[0].get_origin()
    )
    assert shift.elems == pytest.approx((-0.090, -0.168, 0.0), abs=1e-1)


def test_index_after_search(dials_data, run_in_tmp_path):
    """Integrate the beam centre search with the rest of the toolchain

    Do the following:
    1. Take a known good experiment and perturb the beam centre
    2. Run dials.search_beam_centre on the perturbated beam centre and original
    reflection table, check for expected output;
    3. Run dials.index with the found beam centre and check that the expected
    unit cell is obtained and that the RMSDs are smaller than or equal to some
    expected values."""

    insulin = dials_data("insulin_processed")

    # load the original experiment and perturb
    # the beam centre by a small offset
    experiments = load.experiment_list(insulin / "imported.expt", check_format=False)
    original_origin = experiments[0].detector.hierarchy().get_origin()
    shifted_origin = (
        original_origin[0] - 1.3,
        original_origin[1] + 1.5,
        original_origin[2],
    )
    experiments[0].detector.hierarchy().set_local_frame(
        experiments[0].detector.hierarchy().get_fast_axis(),
        experiments[0].detector.hierarchy().get_slow_axis(),
        shifted_origin,
    )
    assert experiments[0].detector.hierarchy().get_origin() == shifted_origin
    experiments.as_file(run_in_tmp_path / "shifted.expt")

    # search the beam centre
    search_beam_position.run(
        [
            str(run_in_tmp_path / "shifted.expt"),
            str(insulin / "strong.refl"),
        ]
    )
    assert run_in_tmp_path.joinpath("optimised.expt").is_file()

    # check we can actually index the resulting optimized experiments
    expected_unit_cell = uctbx.unit_cell(
        (67.655, 67.622, 67.631, 109.4583, 109.4797, 109.485)
    )
    expected_rmsds = (0.3, 0.3, 0.005)
    expected_hall_symbol = " P 1"
    run_indexing(
        insulin / "strong.refl",
        run_in_tmp_path / "optimised.expt",
        run_in_tmp_path,
        [],
        expected_unit_cell,
        expected_rmsds,
        expected_hall_symbol,
    )


def test_search_single(dials_data, run_in_tmp_path):
    """Perform a beam-centre search and check that the output is sane.

    Do the following:
    1. Run dials.search_beam_centre on a single experiment and pickled
    reflection table, as output by dials.find_spots;
      a) Check that the program exits correctly;
      b) Check that it produces the expected output experiment.
    2. Check that the beam centre search has resulted in the expected shift
    in detector origin.
    """

    insulin = dials_data("insulin_processed")
    refl_path = insulin / "strong.refl"
    experiments_path = insulin / "imported.expt"

    search_beam_position.run([str(experiments_path), str(refl_path)])
    assert run_in_tmp_path.joinpath("optimised.expt").is_file()

    experiments = load.experiment_list(experiments_path, check_format=False)
    original_imageset = experiments.imagesets()[0]
    optimized_experiments = load.experiment_list("optimised.expt", check_format=False)
    detector_1 = original_imageset.get_detector()
    detector_2 = optimized_experiments.detectors()[0]
    shift = scitbx.matrix.col(detector_1[0].get_origin()) - scitbx.matrix.col(
        detector_2[0].get_origin()
    )
    assert shift.elems == pytest.approx((-0.165, -0.380, 0.0), abs=1e-1)


def test_search_small_molecule(dials_data, run_in_tmp_path):
    """Perform a beam-centre search on a multi-sequence data set..

    Do the following:
    1. Run dials.search_beam_centre on a single experiment and pickled
    reflection table containing multiple experiment IDs, as output by
    dials.find_spots;
      a) Check that the program exits correctly;
      b) Check that it produces the expected output experiment.
    2. Check that the beam centre search has resulted in the expected shift
    in detector origin.
    """

    data = dials_data("l_cysteine_dials_output")
    experiments_path = data / "imported.expt"
    refl_path = data / "strong.refl"

    search_beam_position.run([os.fspath(experiments_path), os.fspath(refl_path)])
    assert run_in_tmp_path.joinpath("optimised.expt").is_file()

    experiments = load.experiment_list(experiments_path, check_format=False)
    optimised_experiments = load.experiment_list("optimised.expt", check_format=False)
    for old_expt, new_expt in zip(experiments, optimised_experiments):
        # assert that the detector fast/slow axes are unchanged
        # from the input experiments
        # the last experiment actually does have a different detector model
        assert (
            old_expt.detector[0].get_slow_axis() == new_expt.detector[0].get_slow_axis()
        )
        assert (
            old_expt.detector[0].get_fast_axis() == new_expt.detector[0].get_fast_axis()
        )
        shift = scitbx.matrix.col(
            old_expt.detector[0].get_origin()
        ) - scitbx.matrix.col(new_expt.detector[0].get_origin())
        assert shift.elems == pytest.approx((0.091, -1.11, 0), abs=1e-2)


def test_multi_sweep_fixed_rotation(dials_data, run_in_tmp_path):
    data = dials_data("l_cysteine_dials_output")
    experiments_path = data / "imported.expt"
    refl_path = data / "strong.refl"

    search_beam_position.run([str(experiments_path), str(refl_path)])
    assert run_in_tmp_path.joinpath("optimised.expt").is_file()

    experiments = load.experiment_list(experiments_path, check_format=False)
    optimised_experiments = load.experiment_list("optimised.expt", check_format=False)

    for orig_expt, new_expt in zip(experiments, optimised_experiments):
        shift = scitbx.matrix.col(
            orig_expt.detector[0].get_origin()
        ) - scitbx.matrix.col(new_expt.detector[0].get_origin())
        print(shift)
        assert shift.elems == pytest.approx((0.096, -1.111, 0), abs=1e-2)


def test_midpoint_method(tmp_path):
    make_mrc_file(tmp_path)
    imported_file = os.path.join(tmp_path, "imported.expt")

    assert os.path.exists(imported_file), "Error! No imported.expt file."

    dials_cmd = "dials.search_beam_position"
    if os.name == "nt":
        dials_cmd += ".EXE"  # Used for tests on Windows

    cmd = [
        dials_cmd,
        "method=midpoint",
        "per_image=False",
        "exclude_intensity_percent=0.01",
        "intersection_range=0.1,0.9,0.01",
        "midpoint.convolution_width=10",
        "dead_pixel_range_y=223,276",
        "intersection_min_width=10",
        "color_cutoff=20",
        "json=beam_position.json",
        imported_file,
    ]

    result = subprocess.run(cmd, cwd=tmp_path, capture_output=True, text=True)
    check_output(
        result,
        x0=401.57,
        y0=131.00,
        image_index=0,
        imageset_index=0,
        out_fig="beam_position_detector_0.png",
        method="midpoint",
        cwd=tmp_path,
    )

    shutil.rmtree(tmp_path)


@pytest.mark.xdist_group("shared_resource")  # Run within the same process
def test_maximum_method(tmp_path):
    make_mrc_file(tmp_path)
    imported_file = os.path.join(tmp_path, "imported.expt")
    assert os.path.exists(imported_file), "Error! No imported.expt file."

    dials_cmd = "dials.search_beam_position"
    if os.name == "nt":
        dials_cmd += ".EXE"  # Used for tests on Windows

    cmd = [
        dials_cmd,
        "method=maximum",
        "per_image=False",
        "maximum.bad_pixel_threshold=30",
        "maximum.convolution_width=5",
        "bin_width=20",
        "bin_step=10",
        "color_cutoff=20",
        "json=beam_position.json",
        imported_file,
    ]

    result = subprocess.run(cmd, cwd=tmp_path, capture_output=True, text=True)
    check_output(
        result,
        x0=394.0,
        y0=221.0,
        image_index=1,
        imageset_index=0,
        out_fig="beam_position_detector_0.png",
        method="maximum",
        cwd=tmp_path,
    )

    shutil.rmtree(tmp_path)


def test_inversion_method(tmp_path):
    make_mrc_file(tmp_path)

    imported_file = os.path.join(tmp_path, "imported.expt")
    assert os.path.exists(imported_file), "Error! No imported.expt file."

    dials_cmd = "dials.search_beam_position"
    if os.name == "nt":
        dials_cmd += ".EXE"  # Used for tests on Windows

    cmd = [
        dials_cmd,
        "method=inversion",
        "per_image=False",
        "inversion.bad_pixel_threshold=20",
        "color_cutoff=20",
        "json=beam_position.json",
        imported_file,
    ]

    result = subprocess.run(cmd, cwd=tmp_path, capture_output=True, text=True)
    check_output(
        result,
        x0=396.0,
        y0=257.0,
        image_index=2,
        imageset_index=0,
        out_fig="beam_position_detector_0.png",
        method="inversion",
        cwd=tmp_path,
    )

    shutil.rmtree(tmp_path)


def make_mrc_file(path):
    """
    Generates three MRC images to test search_beam_position projection
    methods. The first image is for testing the midpoint method, the second is
    for the maximum method, and the third is for the inverson method
    """

    nx = 700
    ny = 500
    images = []
    beam_positions = [(372, 247), (402, 221), (-3780, -2520)]
    for i in range(3):
        x0, y0 = beam_positions[i]
        images.append(generate_image(nx, ny, x0=x0, y0=y0, seed=i))
    images = np.array(images, dtype=np.float32)

    images[0][225:275, :] = 0  # Add the central stripe to image 1
    images[2][225:275, :] = 0  # Add the central stripe to image 3

    # Add a Friedel pair to image 3
    delta_x = 70
    delta_y = 120

    x0, y0 = 378, 252
    friedel_01 = Gaussian(
        x0=x0 + delta_x, y0=y0 + delta_y, sigma_x=5, sigma_y=5, height=20
    )
    friedel_02 = Gaussian(
        x0=x0 - delta_x, y0=y0 - delta_y, sigma_x=5, sigma_y=5, height=20
    )
    images[2] += gauss(nx, ny, friedel_01)
    images[2] += gauss(nx, ny, friedel_02)

    filename = os.path.join(path, "images.mrc")

    save_to_mrc(images, filename)

    cmd = "dials.import"
    if os.name == "nt":
        cmd += ".EXE"  # Used for tests on Windows

    subprocess.run([cmd, filename], cwd=path, shell=False, text=True)

    return filename

    # os.remove(filename)
    # os.remove(os.join(test_path, "imported.expt"))
    # os.remove(os.join(test_path, "dials.import.log"))


def check_output(
    result,
    x0,
    y0,
    image_index,
    imageset_index,
    out_fig=None,
    method="maximum",
    cwd=None,
):
    out_lines = result.stdout.split("\n")

    assert result.stderr == "", f"Error in subprocess: {result.stderr}"

    # We assume the line in the output where we print the beam position is
    # underneath the detector number (e.g. '[Detector ...]')
    position_line = None
    for index, line in enumerate(out_lines):
        if "[Detector" in line and index < len(out_lines) - 1:
            position_line = out_lines[index + 1]

    msg = f"Function for {method} method prints the output in a wrong format"
    assert position_line is not None, msg

    x, y = [float(i) for i in position_line.strip().split(" ")]
    print(f"Checking for the {method} method")
    print("x, y =", x, y)
    x0_OK = abs(x - x0) < 1
    y0_OK = abs(y - y0) < 1
    msg = f"Wrong beam position from the {method} method \n"
    msg += f"  Got x = {x}, expected {x0} (Position OK? {x0_OK}) \n"
    msg += f"  Got y = {y}, expected {y0} (Position OK? {y0_OK}) \n"

    assert x0_OK and y0_OK, msg

    out_fig = os.path.join(cwd, out_fig)

    assert os.path.exists(out_fig), f"No output figure generated {out_fig}"
    os.remove(out_fig)


def generate_image(nx, ny, x0=500, y0=400, seed=1, num_peaks=100, num_bad_pixels=5):
    xs, ys, widths, heights = get_random_gaussian_data()

    x_pos = xs[seed]
    y_pos = ys[seed]
    hs = heights[seed]
    sigmas = widths[seed]

    image = np.zeros((ny, nx), dtype=np.float64)

    # First, add the central electron beam
    peak = Gaussian(x0=x0, y0=y0, sigma_x=15, sigma_y=15, height=30)
    image += gauss(nx, ny, peak)

    # Add random "diffraction" peaks as narrow Gaussians
    for i in range(len(x_pos)):
        peak = Gaussian(
            x0=x_pos[i], y0=y_pos[i], sigma_x=sigmas[i], sigma_y=sigmas[i], height=hs[i]
        )
        image += gauss(nx, ny, peak)

    return image


def gauss(nx, ny, params):
    y, x = np.meshgrid(np.arange(ny), np.arange(nx), indexing="ij")

    x_rot = -((x - params.x0) ** 2) / (2 * params.sigma_x**2)
    y_rot = -((y - params.y0) ** 2) / (2 * params.sigma_y**2)

    gaussian = params.height * np.exp(x_rot + y_rot)

    return gaussian


def save_to_mrc(images, filename="images.mrc"):
    work_dir = os.path.dirname(filename)
    os.makedirs(work_dir, exist_ok=True)

    with mrcfile.new(filename, overwrite=True) as mrc:
        mrc.set_data(images)

        mrc.voxel_size = 1.0  # example voxel size in angstroms
        mrc.header.nxstart, mrc.header.nystart, mrc.header.nzstart = 0, 0, 0


def get_random_gaussian_data():
    """
    Return data necessary to generate MRC file with diffraction images.
    """

    # This all is random generated data. The reason we keep it is because
    # the generated random numbers can possibly cange with a new Python
    # version.

    # The x positions of the random Gaussians on the first diffraction image
    x1 = [
        394,
        265,
        310,
        223,
        97,
        545,
        317,
        483,
        323,
        565,
        266,
        14,
        684,
        341,
        195,
        145,
        82,
        111,
        560,
        208,
        601,
        394,
        188,
        627,
        91,
        39,
        699,
        282,
        695,
        593,
        365,
        498,
        248,
        225,
        436,
        224,
        547,
        127,
        122,
        37,
        126,
        695,
        103,
        73,
        184,
        610,
        204,
        583,
        59,
        350,
        452,
        697,
        318,
        157,
        80,
        287,
        360,
        604,
        135,
        666,
        342,
        458,
        424,
        428,
        169,
        457,
        618,
        506,
        51,
        192,
        85,
        411,
        3,
        691,
        121,
        286,
        487,
        83,
        118,
        533,
        117,
        43,
        571,
        43,
        671,
        655,
        182,
        9,
        345,
        346,
        276,
        296,
        132,
        49,
        72,
        306,
        492,
        127,
        436,
        170,
        89,
        226,
        403,
        500,
        433,
    ]

    # The x positions of the random Gaussians on the second diffraction image
    x2 = [
        137,
        64,
        460,
        214,
        399,
        2,
        605,
        31,
        9,
        432,
        448,
        353,
        296,
        657,
        123,
        432,
        310,
        517,
        491,
        424,
        690,
        679,
        533,
        30,
        607,
        172,
        204,
        237,
        361,
        623,
        524,
        210,
        373,
        423,
        354,
        627,
        235,
        185,
        261,
        85,
        287,
        189,
        163,
        672,
        465,
        116,
        431,
        259,
        620,
        230,
        164,
        557,
        528,
        31,
        675,
        128,
        313,
        305,
        258,
        604,
        471,
        521,
        101,
        605,
        303,
        626,
        18,
        576,
        218,
        560,
        496,
        66,
        173,
        274,
        261,
        298,
        500,
        328,
        150,
        629,
        78,
        83,
        577,
        283,
        12,
        423,
        192,
        431,
        697,
        445,
        301,
        322,
        40,
        610,
        408,
        615,
        220,
        480,
        554,
        369,
        458,
        347,
        42,
        310,
        626,
    ]

    # The x positions of the random Gaussians on the third diffraction image
    x3 = [
        57,
        173,
        315,
        36,
        653,
        521,
        514,
        372,
        433,
        573,
        180,
        522,
        186,
        537,
        607,
        165,
        670,
        510,
        677,
        359,
        467,
        170,
        491,
        575,
        602,
        500,
        77,
        195,
        668,
        698,
        139,
        61,
        58,
        688,
        25,
        382,
        188,
        603,
        155,
        630,
        292,
        459,
        270,
        157,
        676,
        24,
        530,
        527,
        265,
        668,
        686,
        165,
        237,
        252,
        256,
        638,
        285,
        36,
        524,
        102,
        107,
        464,
        217,
        444,
        605,
        595,
        374,
        625,
        315,
        103,
        462,
        474,
        75,
        313,
        502,
        401,
        404,
        125,
        656,
        25,
        509,
        144,
        538,
        503,
        237,
        561,
        599,
        72,
        150,
        88,
        38,
        237,
        449,
        159,
        98,
        263,
        240,
        587,
        30,
        402,
        414,
        663,
        338,
        230,
        682,
    ]

    # The y positions of the random Gaussians on the first diffraction image
    y1 = [
        388,
        494,
        495,
        258,
        316,
        361,
        50,
        286,
        312,
        244,
        31,
        47,
        320,
        124,
        469,
        411,
        163,
        154,
        170,
        493,
        147,
        162,
        96,
        336,
        347,
        431,
        200,
        267,
        301,
        140,
        42,
        300,
        8,
        190,
        417,
        23,
        308,
        325,
        200,
        310,
        245,
        11,
        427,
        331,
        31,
        51,
        133,
        86,
        403,
        271,
        340,
        209,
        332,
        287,
        171,
        69,
        312,
        484,
        366,
        41,
        81,
        193,
        16,
        395,
        228,
        270,
        386,
        166,
        414,
        280,
        428,
        486,
        109,
        270,
        311,
        352,
        437,
        11,
        441,
        418,
        446,
        20,
        161,
        433,
        253,
        222,
        106,
        70,
        404,
        399,
        83,
        184,
        150,
        157,
        154,
        212,
        242,
        245,
        19,
        320,
        33,
        31,
        284,
        403,
        42,
    ]

    # The y positions of the random Gaussians on the second diffraction image
    y2 = [
        291,
        130,
        241,
        48,
        221,
        356,
        483,
        11,
        480,
        371,
        480,
        118,
        474,
        51,
        380,
        259,
        145,
        201,
        124,
        340,
        377,
        260,
        430,
        240,
        296,
        257,
        276,
        207,
        235,
        490,
        414,
        218,
        291,
        248,
        0,
        169,
        325,
        440,
        16,
        444,
        127,
        176,
        130,
        139,
        359,
        12,
        407,
        460,
        221,
        9,
        228,
        426,
        230,
        202,
        323,
        495,
        36,
        380,
        66,
        419,
        87,
        19,
        105,
        99,
        258,
        446,
        80,
        400,
        136,
        176,
        393,
        371,
        85,
        388,
        188,
        120,
        69,
        20,
        424,
        300,
        292,
        487,
        273,
        55,
        314,
        58,
        122,
        82,
        123,
        466,
        281,
        51,
        13,
        163,
        32,
        496,
        401,
        338,
        106,
        41,
        46,
        481,
        167,
        125,
        296,
    ]

    # The y positions of the random Gaussians on the third diffraction image
    y3 = [
        46,
        376,
        128,
        297,
        201,
        486,
        137,
        238,
        456,
        90,
        166,
        184,
        457,
        464,
        181,
        488,
        271,
        256,
        452,
        290,
        249,
        448,
        158,
        265,
        208,
        262,
        401,
        381,
        25,
        448,
        437,
        216,
        185,
        12,
        20,
        130,
        267,
        22,
        497,
        321,
        172,
        282,
        386,
        242,
        351,
        229,
        299,
        167,
        134,
        9,
        29,
        87,
        260,
        119,
        41,
        363,
        269,
        196,
        370,
        10,
        111,
        158,
        350,
        217,
        26,
        92,
        9,
        187,
        9,
        53,
        30,
        106,
        2,
        478,
        98,
        366,
        454,
        41,
        200,
        316,
        148,
        408,
        443,
        428,
        80,
        218,
        372,
        182,
        412,
        466,
        65,
        362,
        88,
        276,
        168,
        86,
        206,
        373,
        414,
        261,
        129,
        241,
        365,
        272,
        195,
    ]

    # The widths of the random Gaussians on the first diffraction image
    w1 = [
        0.90121953,
        0.56014725,
        0.52893726,
        0.22534633,
        0.81946232,
        0.82919551,
        0.75684857,
        0.19063109,
        0.67638106,
        0.49842402,
        0.82452504,
        0.74773422,
        0.10102854,
        0.75725072,
        0.61075967,
        0.58869096,
        0.88757856,
        0.59614052,
        0.83302018,
        0.8194844,
        0.50049012,
        0.61808591,
        0.83932049,
        0.33405059,
        0.78184353,
        0.17222349,
        0.85407871,
        0.83046584,
        0.84273697,
        0.5055068,
        0.39189236,
        0.66713261,
        0.7583403,
        0.81512467,
        0.15597321,
        0.83540827,
        0.71253022,
        0.26966856,
        0.18237867,
        0.11947286,
        0.28951886,
        0.58981187,
        0.3339279,
        0.37095538,
        0.55327285,
        0.72939081,
        0.42269977,
        0.72792736,
        0.70858161,
        0.32561861,
        0.25733754,
        0.91008825,
        0.42143651,
        0.72174048,
        0.7651145,
        0.31611737,
        0.35909012,
        0.90274919,
        0.37926782,
        0.1013718,
        0.31547436,
        0.73931498,
        0.46204037,
        0.6959763,
        0.15749768,
        0.89645096,
        0.10006225,
        0.38085242,
        0.84044793,
        0.99947832,
        0.75299892,
        0.7106635,
        0.11286218,
        0.65085577,
        0.68433888,
        0.98066367,
        0.93049276,
        0.34724246,
        0.33083989,
        0.68570491,
        0.23899612,
        0.28516217,
        0.9518523,
        0.77422134,
        0.74112189,
        0.43520434,
        0.43802103,
        0.23592086,
        0.43048133,
        0.6584843,
        0.23447443,
        0.45531278,
        0.20340677,
        0.26162816,
        0.46285223,
        0.19778627,
        0.40337249,
        0.20440726,
        0.3717128,
        0.60803865,
        0.82699531,
        0.44633071,
        0.56710336,
        0.62639983,
        0.43149831,
    ]

    # The widths of the random Gaussians on the second diffraction image
    w2 = [
        0.86269036,
        0.20612683,
        0.68643368,
        0.53907099,
        0.6466942,
        0.50084847,
        0.19200444,
        0.12290127,
        0.89311047,
        0.12613671,
        0.54623102,
        0.70916369,
        0.11934073,
        0.26731564,
        0.39942567,
        0.84703212,
        0.62882255,
        0.63010203,
        0.76935759,
        0.25570666,
        0.43723272,
        0.19713814,
        0.45392959,
        0.13913856,
        0.45423972,
        0.30424361,
        0.92807779,
        0.5623945,
        0.91877951,
        0.75646016,
        0.21632448,
        0.95474734,
        0.59896107,
        0.83201636,
        0.58463092,
        0.51233212,
        0.25949013,
        0.18244127,
        0.85757035,
        0.11502157,
        0.34177756,
        0.36125589,
        0.57464236,
        0.68338252,
        0.38980159,
        0.38080443,
        0.26923537,
        0.75713701,
        0.8353357,
        0.45758383,
        0.73414507,
        0.29853977,
        0.3008762,
        0.70741232,
        0.48366944,
        0.29091852,
        0.87266297,
        0.24238665,
        0.10763224,
        0.29583131,
        0.84512184,
        0.44017195,
        0.61604443,
        0.54311568,
        0.54981993,
        0.46209224,
        0.28076771,
        0.22162883,
        0.70700956,
        0.92300018,
        0.98695741,
        0.13636041,
        0.91988901,
        0.39901013,
        0.40496642,
        0.88067784,
        0.62196135,
        0.46593903,
        0.21251553,
        0.80371324,
        0.59524728,
        0.34005128,
        0.93256048,
        0.8083132,
        0.70337048,
        0.84350799,
        0.80696092,
        0.20400236,
        0.24305029,
        0.96822442,
        0.32806912,
        0.28685316,
        0.10945554,
        0.50486008,
        0.15777278,
        0.51027219,
        0.65602619,
        0.42024428,
        0.37658901,
        0.83797934,
        0.68680948,
        0.30468361,
        0.26813157,
        0.40090007,
        0.82681093,
    ]

    # The widths of the random Gaussians on the third diffraction image
    w3 = [
        0.1763848,
        0.82806584,
        0.64534975,
        0.71312157,
        0.8233025,
        0.43485427,
        0.91108702,
        0.9386619,
        0.89565811,
        0.31251107,
        0.25622988,
        0.99881521,
        0.50107385,
        0.91763027,
        0.42569774,
        0.77860366,
        0.32491159,
        0.56385298,
        0.50921991,
        0.75329668,
        0.69299677,
        0.91714645,
        0.37296707,
        0.55662102,
        0.38065989,
        0.42993842,
        0.83836,
        0.19554123,
        0.34578976,
        0.92425687,
        0.33925708,
        0.90817109,
        0.42419663,
        0.17461286,
        0.75648116,
        0.21500842,
        0.72234588,
        0.81474213,
        0.13264209,
        0.76821366,
        0.53988492,
        0.78937294,
        0.46166859,
        0.96290564,
        0.38467533,
        0.8095684,
        0.80295366,
        0.22944652,
        0.64547015,
        0.72952643,
        0.32767023,
        0.18628475,
        0.92475817,
        0.74248073,
        0.63211361,
        0.42384612,
        0.77561001,
        0.46785378,
        0.17903922,
        0.26358407,
        0.1219961,
        0.58198001,
        0.91590641,
        0.56035121,
        0.89329334,
        0.9285547,
        0.56728872,
        0.36062269,
        0.88691231,
        0.37533012,
        0.46960694,
        0.90058772,
        0.35592342,
        0.7512085,
        0.20419153,
        0.51690194,
        0.20955653,
        0.65488732,
        0.96370171,
        0.6932437,
        0.42182616,
        0.43742125,
        0.52980359,
        0.71220122,
        0.53991053,
        0.72605668,
        0.83606793,
        0.25852095,
        0.47509124,
        0.90595032,
        0.9765232,
        0.70365909,
        0.57159529,
        0.98801767,
        0.56485123,
        0.90540698,
        0.88702139,
        0.23028852,
        0.63564776,
        0.1482151,
        0.73941407,
        0.42409609,
        0.77105148,
        0.65912413,
        0.97951066,
    ]

    # The heights of the random Gaussians on the first diffraction image
    h1 = [
        71,
        444,
        627,
        173,
        961,
        4,
        10,
        13,
        6,
        16,
        17,
        12,
        15,
        2,
        7,
        2,
        15,
        3,
        17,
        17,
        19,
        9,
        1,
        2,
        4,
        17,
        16,
        6,
        13,
        20,
        3,
        6,
        3,
        10,
        4,
        20,
        0,
        18,
        3,
        5,
        1,
        19,
        7,
        13,
        1,
        12,
        15,
        6,
        5,
        19,
        15,
        16,
        8,
        14,
        17,
        15,
        11,
        19,
        13,
        6,
        20,
        18,
        18,
        1,
        5,
        17,
        1,
        14,
        13,
        20,
        0,
        10,
        0,
        6,
        9,
        3,
        20,
        14,
        20,
        11,
        0,
        8,
        18,
        19,
        14,
        17,
        9,
        10,
        2,
        1,
        18,
        17,
        7,
        16,
        10,
        17,
        10,
        15,
        4,
        20,
        6,
        3,
        14,
        6,
        8,
    ]

    # The heights of the random Gaussians on the second diffraction image
    h2 = [
        812,
        809,
        837,
        944,
        815,
        7,
        10,
        17,
        6,
        7,
        7,
        14,
        17,
        9,
        16,
        6,
        15,
        1,
        12,
        17,
        14,
        5,
        15,
        19,
        5,
        0,
        17,
        18,
        17,
        12,
        17,
        15,
        16,
        13,
        19,
        0,
        18,
        17,
        2,
        0,
        19,
        5,
        5,
        9,
        15,
        10,
        3,
        6,
        0,
        1,
        13,
        20,
        20,
        10,
        9,
        1,
        9,
        18,
        1,
        18,
        19,
        11,
        13,
        12,
        10,
        9,
        10,
        13,
        12,
        17,
        7,
        4,
        6,
        16,
        3,
        19,
        3,
        12,
        3,
        12,
        18,
        9,
        14,
        9,
        2,
        1,
        18,
        5,
        3,
        17,
        15,
        10,
        9,
        10,
        10,
        8,
        17,
        5,
        7,
        2,
        20,
        9,
        18,
        17,
        2,
    ]

    # The heights of the random Gaussians on the third diffraction image
    h3 = [
        885,
        904,
        651,
        471,
        910,
        14,
        0,
        12,
        5,
        0,
        16,
        17,
        13,
        11,
        14,
        14,
        8,
        11,
        14,
        17,
        10,
        8,
        16,
        19,
        6,
        19,
        0,
        18,
        7,
        16,
        6,
        1,
        7,
        2,
        0,
        5,
        12,
        7,
        11,
        3,
        9,
        1,
        19,
        2,
        3,
        4,
        15,
        10,
        13,
        4,
        4,
        20,
        1,
        2,
        19,
        13,
        4,
        3,
        3,
        7,
        14,
        12,
        6,
        18,
        16,
        15,
        3,
        11,
        13,
        0,
        15,
        19,
        11,
        7,
        11,
        11,
        3,
        10,
        3,
        1,
        14,
        15,
        13,
        12,
        8,
        2,
        3,
        17,
        2,
        20,
        12,
        10,
        3,
        13,
        16,
        14,
        11,
        14,
        5,
        8,
        13,
        17,
        2,
        12,
        20,
    ]

    xs = [x1, x2, x3]  # X positions of the random Gaussians
    ys = [y1, y2, y3]  # Y positions of the random Gaussians
    widths = [w1, w2, w3]  # Widths of the random Gaussians
    heights = [h1, h2, h3]  # Heights of the random Gaussians

    return xs, ys, widths, heights
