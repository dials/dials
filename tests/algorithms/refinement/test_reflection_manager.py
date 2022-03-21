from __future__ import annotations

import math

import pytest

from dxtbx.model.experiment_list import ExperimentListFactory

from dials.algorithms.refinement.reflection_manager import ReflectionManager
from dials.array_family import flex


def test_scan_margin(dials_data):

    # Use 4 scan data for this test
    data_dir = dials_data("l_cysteine_dials_output", pathlib=True)
    experiments = ExperimentListFactory.from_json_file(
        data_dir / "indexed.expt", check_format=False
    )
    reflections = flex.reflection_table.from_file(data_dir / "indexed.refl")
    orig_phi = reflections["xyzobs.mm.value"].parts()[2]

    # Reflection Manager works on predictions, but this dataset has none, so
    # need to set the predictions flag
    reflections.set_flags(
        flex.bool(len(reflections), True), reflections.flags.predicted
    )

    # Create a reflection manager without trimming scan margins
    refman = ReflectionManager(reflections, experiments)
    refman.finalise()
    refs1 = refman.get_matches()
    phi1 = refs1["xyzobs.mm.value"].parts()[2]

    # Create a reflection manager with 1 degree width scan margins
    margin = 1.0
    refman2 = ReflectionManager(reflections, experiments, scan_margin=margin)
    refman2.finalise()
    refs2 = refman2.get_matches()
    phi2 = refs2["xyzobs.mm.value"].parts()[2]

    # Check zero scan margins do not trim
    assert min(orig_phi) == min(phi1)
    assert max(orig_phi) == max(phi1)

    # Check 1 degree scan margin trims approximately 1 degree
    assert min(phi2) == pytest.approx(min(phi1) + math.radians(margin), abs=1e-3)
    assert max(phi2) == pytest.approx(max(phi1) - math.radians(margin), abs=1e-3)
