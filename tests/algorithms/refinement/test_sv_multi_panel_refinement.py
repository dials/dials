from __future__ import annotations

import shutil
import subprocess

from dials.algorithms.refinement.engine import Journal


def test_scan_varying_refinement_of_a_multiple_panel_detector(dials_data, tmp_path):
    from dials.array_family import flex

    data_dir = dials_data("refinement_test_data", pathlib=True)

    result = subprocess.run(
        [
            shutil.which("dials.refine"),
            data_dir / "I23_24_panel.json",
            data_dir / "I23_24_panel.pickle",
            "scan_varying=true",
            "history=history.json",
            "outlier.separate_blocks=False",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    # there are plenty of things we could do with the refinement history, but
    # here just check that final RMSDs are low enough
    history = Journal.from_json_file(tmp_path / "history.json")
    final_rmsd = history["rmsd"][-1]
    assert final_rmsd[0] < 0.05
    assert final_rmsd[1] < 0.04
    assert final_rmsd[2] < 0.0002

    # also check that the used_in_refinement flag got set correctly
    rt = flex.reflection_table.from_file(tmp_path / "refined.refl")
    uir = rt.get_flags(rt.flags.used_in_refinement)
    assert uir.count(True) == history["num_reflections"][-1]
