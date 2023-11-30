"""
Test command line program dials.refine by running a job with saved data and
comparing with expected output.

This serves as a high level test that not only checks whether refinement works,
but also that the command line program is functioning and that the output models
have not changed format and so on.
"""


from __future__ import annotations

import math
import random
import shutil
import subprocess

import numpy as np
import pytest
from annlib_ext import AnnAdaptor

from dxtbx.model import Beam, Crystal, Detector, Experiment, Goniometer
from dxtbx.model.experiment_list import ExperimentList, ExperimentListFactory

from dials.algorithms.refinement.engine import Journal
from dials.array_family import flex
from dials.command_line.refine import _find_disjoint_sets


def test_i04_weak_data(dials_data, tmp_path):
    data_dir = dials_data("refinement_test_data", pathlib=True)
    experiments_path = data_dir / "i04-weak.json"
    pickle_path = data_dir / "i04-weak.pickle"

    # set some old defaults
    cmd = (
        shutil.which("dials.refine"),
        "close_to_spindle_cutoff=0.05",
        "reflections_per_degree=100",
        "outlier.separate_blocks=False",
        "scan_varying=False",
        "reflections.outlier.nproc=1",
        experiments_path,
        pickle_path,
    )
    result = subprocess.run(cmd, cwd=tmp_path, capture_output=True)
    assert not result.returncode and not result.stderr
    # load results
    reg_exp = ExperimentListFactory.from_json_file(
        data_dir / "i04-weak-regression.json", check_format=False
    )[0]
    ref_exp = ExperimentListFactory.from_json_file(
        tmp_path / "refined.expt", check_format=False
    )[0]

    # test refined models against expected
    assert reg_exp.crystal == ref_exp.crystal
    assert reg_exp.detector == ref_exp.detector
    assert reg_exp.beam == ref_exp.beam

    # test cell parameter esds
    assert ref_exp.crystal.get_cell_parameter_sd() == pytest.approx(
        (0.0009903, 0.0009903, 0.0021227, 0.0, 0.0, 0.0), abs=1e-6
    )
    assert ref_exp.crystal.get_cell_volume_sd() == pytest.approx(23.8051262, abs=1e-6)


def test_scan_varying_refinement_rmsds(dials_data, tmp_path):
    """Run scan-varying refinement, comparing RMSD table with expected values."""

    data_dir = dials_data("refinement_test_data", pathlib=True)
    experiments_path = data_dir / "from-xds.json"
    pickle_path = data_dir / "from-xds-all.pickle"

    # scan-static refinement first to get refined.expt as start point
    result = subprocess.run(
        (
            shutil.which("dials.refine"),
            experiments_path,
            pickle_path,
            "scan_varying=False",
            "reflections_per_degree=50",
            "outlier.algorithm=null",
            "close_to_spindle_cutoff=0.05",
        ),
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    # NB Requesting corrgram.pdf and history.json exercises
    # https://github.com/dials/dials/issues/1923
    result = subprocess.run(
        (
            shutil.which("dials.refine"),
            "refined.expt",
            pickle_path,
            "scan_varying=true",
            "output.history=history.json",
            "correlation_plot.filename=corrgram.pdf",
            "reflections_per_degree=50",
            "outlier.algorithm=null",
            "close_to_spindle_cutoff=0.05",
            "crystal.orientation.smoother.interval_width_degrees=36.0",
            "crystal.unit_cell.smoother.interval_width_degrees=36.0",
            "set_scan_varying_errors=True",
            "reflections.outlier.nproc=1",
        ),
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    # load and check results
    history = Journal.from_json_file(tmp_path / "history.json")

    expected_rmsds = [
        (0.088488398, 0.114583571, 0.001460382),
        (0.080489334, 0.086406517, 0.001284069),
        (0.078835086, 0.086052630, 0.001195882),
        (0.077476911, 0.086194611, 0.001161143),
        (0.076755840, 0.086090630, 0.001157239),
        (0.076586376, 0.085939462, 0.001155641),
        (0.076603722, 0.085878953, 0.001155065),
        (0.076611382, 0.085862959, 0.001154863),
        (0.076608732, 0.085856935, 0.001154384),
        (0.076605731, 0.085852271, 0.001153858),
        (0.076604576, 0.085852318, 0.001153643),
        (0.076603981, 0.085854175, 0.001153594),
    ]
    for a, b in zip(history["rmsd"], expected_rmsds):
        assert a == pytest.approx(b, abs=1e-6)

    # check that the used_in_refinement flag got set correctly
    rt = flex.reflection_table.from_file(tmp_path / "refined.refl")
    uir = rt.get_flags(rt.flags.used_in_refinement)
    assert uir.count(True) == history["num_reflections"][-1]


def test_scan_varying_with_automated_outlier_rejection_block_width_interval_width(
    dials_data, tmp_path
):
    """Strict check for scan-varying refinement using automated outlier rejection
    block width and interval width setting"""

    data_dir = dials_data("refinement_test_data", pathlib=True)
    experiments_path = data_dir / "from-xds.json"
    pickle_path = data_dir / "from-xds-all.pickle"

    result = subprocess.run(
        (
            shutil.which("dials.refine"),
            experiments_path,
            pickle_path,
            "scan_varying=true",
            "max_iterations=5",
            "output.history=history.json",
            "crystal.orientation.smoother.interval_width_degrees=auto",
            "crystal.unit_cell.smoother.interval_width_degrees=auto",
            "reflections.outlier.nproc=1",
        ),
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    # load and check results
    history = Journal.from_json_file(tmp_path / "history.json")

    expected_rmsds = [
        [0.619507829, 0.351326044, 0.006955399],
        [0.174024575, 0.113486044, 0.004704006],
        [0.098351363, 0.084052519, 0.002660408],
        [0.069202909, 0.072796782, 0.001451734],
        [0.064305277, 0.071560831, 0.001165639],
        [0.062955462, 0.071315612, 0.001074453],
    ]
    for a, b in zip(history["rmsd"], expected_rmsds):
        assert a == pytest.approx(b, abs=1e-6)

    # check the refined unit cell
    ref_exp = ExperimentListFactory.from_json_file(
        tmp_path / "refined.expt", check_format=False
    )[0]
    unit_cell = ref_exp.crystal.get_unit_cell().parameters()
    assert unit_cell == pytest.approx(
        [42.27482, 42.27482, 39.66893, 90.00000, 90.00000, 90.00000], abs=1e-3
    )

    refined_refl = flex.reflection_table.from_file(tmp_path / "refined.refl")
    # re-predict reflections using the refined experiments
    predicted = flex.reflection_table.from_predictions_multi([ref_exp])

    matched, reference, unmatched = predicted.match_with_reference(refined_refl)
    # assert most refined reflections are matched with predictions
    assert reference.size() > (0.997 * refined_refl.size())

    # second check with nearest neighbour matching that the predictions match up
    ann = AnnAdaptor(data=predicted["xyzcal.px"].as_double(), dim=3, k=1)
    ann.query(refined_refl["xyzcal.px"].as_double())
    assert (ann.distances < 0.5).count(True) > (0.998 * refined_refl.size())


@pytest.mark.parametrize("fix", ["cell", "orientation"])
def test_scan_varying_with_fixed_crystal(fix, dials_data, tmp_path):
    location = dials_data("multi_crystal_proteinase_k", pathlib=True)
    refls = location / "reflections_1.pickle"
    expts = location / "experiments_1.json"

    result = subprocess.run(
        [shutil.which("dials.refine"), expts, refls, f"crystal.fix={fix}"],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr


def test_scan_varying_missing_segments_multi_crystal(dials_data, tmp_path):
    # https://github.com/dials/dials/issues/1053
    location = dials_data("i19_1_pdteet_index", pathlib=True)
    refls = location / "indexed.refl"
    expts = location / "indexed.expt"

    # first remove some reflections to keep dials.refine sharp
    data = flex.reflection_table.from_file(refls)

    # prune only indexed reflections

    data = data.select(data.get_flags(data.flags.indexed))
    z = data["xyzobs.px.value"].parts()[2]

    # take overlapping subsets of reflections from two lattices which do not
    # fill the entire scan

    sel = ((data["id"] == 0) & (z > 200)) | ((data["id"] == 1) & (z < 650))
    trimmed = data.select(sel)

    trimmed_filename = tmp_path / "indexed_trim.refl"
    trimmed.as_file(trimmed_filename)

    result = subprocess.run(
        [shutil.which("dials.refine"), expts, trimmed_filename],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    el = ExperimentListFactory.from_json_file(
        tmp_path / "refined.expt", check_format=False
    )

    # verify scans start and finish correctly
    image_ranges = [e.scan.get_image_range() for e in el]

    assert image_ranges[0][1] == 850
    assert image_ranges[1][0] == 1


@pytest.mark.parametrize("gcb", ["None", "3000"])
def test_scan_varying_multi_scan_one_crystal(gcb, dials_data, tmp_path):
    # https://github.com/dials/dials/issues/994
    location = dials_data("l_cysteine_dials_output", pathlib=True)
    refls = location / "indexed.refl"
    expts = location / "indexed.expt"

    # Set options for quick rather than careful refinement
    result = subprocess.run(
        [
            shutil.which("dials.refine"),
            expts,
            refls,
            "output.history=history.json",
            "outlier.algorithm=tukey",
            "max_iterations=3",
            "unit_cell.smoother.interval_width_degrees=56",
            "orientation.smoother.interval_width_degrees=56",
            "gradient_calculation_blocksize=" + gcb,
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    el = ExperimentListFactory.from_json_file(
        tmp_path / "refined.expt", check_format=False
    )

    # Crystal has been copied into each experiment for scan-varying refinement
    assert len(el.crystals()) == 4

    # load and check results
    history = Journal.from_json_file(tmp_path / "history.json")

    expected_rmsds = [
        (0.1401658782847504, 0.2225931584837884, 0.002349912655443814),
        (0.12060230585178289, 0.1585977879739876, 0.002114318828411418),
        (0.10970832317567975, 0.1348574975434352, 0.001955034565537597),
        (0.10373159352273859, 0.12827852889951505, 0.0017901404193256304),
    ]
    for a, b in zip(history["rmsd"], expected_rmsds):
        assert a == pytest.approx(b, abs=1e-6)


def test_find_disjoint_sets():
    # Create some experiments with a shared Beam
    beam = Beam()
    shared_beam = []
    for i in range(random.randint(1, 5)):
        crystal = Crystal((10, 0, 0, 0, 10, 0, 0, 0, 10), "P1")
        goniometer = Goniometer()
        detector = Detector()
        shared_beam.append(
            Experiment(
                beam=beam, goniometer=goniometer, crystal=crystal, detector=detector
            )
        )

    # Create some experiments with a shared Crystal
    crystal = Crystal((10, 0, 0, 0, 10, 0, 0, 0, 10), "P1")
    shared_crystal = []
    for i in range(random.randint(1, 5)):
        beam = Beam()
        goniometer = Goniometer()
        detector = Detector()
        shared_crystal.append(
            Experiment(
                beam=beam, goniometer=goniometer, crystal=crystal, detector=detector
            )
        )

    # Create some experiments with a shared Detector
    detector = Detector()
    shared_detector = []
    for i in range(random.randint(1, 5)):
        beam = Beam()
        crystal = Crystal((10, 0, 0, 0, 10, 0, 0, 0, 10), "P1")
        goniometer = Goniometer()
        shared_detector.append(
            Experiment(
                beam=beam, goniometer=goniometer, crystal=crystal, detector=detector
            )
        )

    # Create some experiments with a shared Goniometer
    goniometer = Goniometer()
    shared_goniometer = []
    for i in range(random.randint(1, 5)):
        beam = Beam()
        crystal = Crystal((10, 0, 0, 0, 10, 0, 0, 0, 10), "P1")
        detector = Detector()
        shared_goniometer.append(
            Experiment(
                beam=beam, goniometer=goniometer, crystal=crystal, detector=detector
            )
        )

    # Randomly re-order and make an experiment list
    experiments = shared_beam + shared_crystal + shared_detector + shared_goniometer
    random.shuffle(experiments)
    experiments = ExperimentList(experiments)

    disjoint_sets = _find_disjoint_sets(experiments)

    # We know there are 4 groups
    assert len(disjoint_sets) == 4

    # The length of each group must match one of the known group lengths
    check = {
        len(shared_beam),
        len(shared_crystal),
        len(shared_detector),
        len(shared_goniometer),
    }
    for group in disjoint_sets:
        assert len(group) in check


def test_refinement_of_disjoint_sets(dials_data, tmp_path):
    # Take 22 input experiments that are in 7 disjoint groups
    location = dials_data("polyhedra_narrow_wedges", pathlib=True)
    expts = [
        location / f"sweep_00{i}_experiments.json"
        for i in ["2", "3", "4", "5", "6", "7", "9"]
    ]
    refls = [
        location / f"sweep_00{i}_reflections.pickle"
        for i in ["2", "3", "4", "5", "6", "7", "9"]
    ]
    result = subprocess.run(
        [
            shutil.which("dials.combine_experiments"),
        ]
        + expts
        + refls,
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    # Refine treating the disjoint groups separately and calculate final RMSDs
    result = subprocess.run(
        [
            shutil.which("dials.refine"),
            "combined.expt",
            "combined.refl",
            "scan_varying=False",
            "output.reflections=separate.refl",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    separate_refl = flex.reflection_table.from_file(tmp_path / "separate.refl")
    separate_refl = separate_refl.select(
        separate_refl.get_flags(separate_refl.flags.used_in_refinement)
    )
    x_obs, y_obs, z_obs = separate_refl["xyzobs.mm.value"].parts()
    x_calc, y_calc, z_calc = separate_refl["xyzcal.mm"].parts()
    x_resid2 = (x_calc - x_obs) ** 2
    y_resid2 = (y_calc - y_obs) ** 2
    z_resid2 = (z_calc - z_obs) ** 2
    separate_rmsd_x = []
    separate_rmsd_y = []
    separate_rmsd_z = []
    for i in range(22):
        sel = separate_refl["id"] == i
        n = sel.count(True)
        separate_rmsd_x.append(math.sqrt(flex.sum(x_resid2.select(sel)) / n))
        separate_rmsd_y.append(math.sqrt(flex.sum(y_resid2.select(sel)) / n))
        separate_rmsd_z.append(math.sqrt(flex.sum(z_resid2.select(sel)) / n))

    # Refine all together in one joint refinement and calculate final RMSDs
    result = subprocess.run(
        [
            shutil.which("dials.refine"),
            "combined.expt",
            "combined.refl",
            "scan_varying=False",
            "separate_independent_sets=False",
            "output.reflections=joint.refl",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    joint_refl = flex.reflection_table.from_file(tmp_path / "joint.refl")
    joint_refl = joint_refl.select(
        joint_refl.get_flags(joint_refl.flags.used_in_refinement)
    )
    x_obs, y_obs, z_obs = joint_refl["xyzobs.mm.value"].parts()
    x_calc, y_calc, z_calc = joint_refl["xyzcal.mm"].parts()
    x_resid2 = (x_calc - x_obs) ** 2
    y_resid2 = (y_calc - y_obs) ** 2
    z_resid2 = (z_calc - z_obs) ** 2
    joint_rmsd_x = []
    joint_rmsd_y = []
    joint_rmsd_z = []
    for i in range(22):
        sel = joint_refl["id"] == i
        n = sel.count(True)
        joint_rmsd_x.append(math.sqrt(flex.sum(x_resid2.select(sel)) / n))
        joint_rmsd_y.append(math.sqrt(flex.sum(y_resid2.select(sel)) / n))
        joint_rmsd_z.append(math.sqrt(flex.sum(z_resid2.select(sel)) / n))

    # Check difference in RMSDs within 1/100th of a pixel and 1/100th of a degree.
    # We don't expect them to be exactly the same (https://github.com/dials/dials/pull/2336#issue-1572772140)
    # but we want the two ways to be of equivalent quality.
    assert (
        flex.max(flex.abs(flex.double(separate_rmsd_x) - flex.double(joint_rmsd_x)))
        < 0.01
    )
    assert (
        flex.max(flex.abs(flex.double(separate_rmsd_y) - flex.double(joint_rmsd_y)))
        < 0.01
    )
    assert (
        np.degrees(
            flex.max(flex.abs(flex.double(separate_rmsd_z) - flex.double(joint_rmsd_z)))
        )
        < 0.01
    )
