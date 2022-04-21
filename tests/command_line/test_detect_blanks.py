from __future__ import annotations

import json

import procrunner
import pytest

from dials.array_family import flex
from dials.command_line import detect_blanks


def test_strong(dials_data, capsys, run_in_tmp_path):
    expts_file = dials_data("insulin_processed", pathlib=True) / "imported.expt"
    refl_file = dials_data("insulin_processed", pathlib=True) / "strong.refl"
    refl = flex.reflection_table.from_file(refl_file)
    z = refl["xyzobs.px.value"].parts()[2]
    refl_subset = refl.select((z < 10) | (z > 20))
    refl_subset.as_file(run_in_tmp_path / "mod.refl")

    detect_blanks.run([str(expts_file), str(run_in_tmp_path / "mod.refl")])

    captured = capsys.readouterr()
    assert "Potential blank images: 11 -> 21" in captured.out
    json_file = run_in_tmp_path / "blanks.json"
    assert json_file.is_file()
    with json_file.open() as fh:
        d = json.load(fh)
    assert d["strong"]["blank_regions"] == [[10, 21]]
    assert not d["indexed"]
    assert not d["integrated"]


def test_integrated(dials_data, capsys, run_in_tmp_path):
    expts_file = dials_data("insulin_processed", pathlib=True) / "integrated.expt"
    refl_file = dials_data("insulin_processed", pathlib=True) / "integrated.refl"
    refl = flex.reflection_table.from_file(refl_file)
    z = refl["xyzobs.px.value"].parts()[2]
    sel = z < 40
    refl_subset = refl.select(sel)
    refl_subset["intensity.prf.value"].set_selected(
        z.select(sel) < 10, refl_subset["intensity.prf.value"] * 0.05
    )
    refl_subset.as_file(run_in_tmp_path / "mod.refl")

    detect_blanks.run(
        [
            str(expts_file),
            str(run_in_tmp_path / "mod.refl"),
            "output.json=detect_blanks.json",
        ]
    )

    captured = capsys.readouterr()
    assert (
        """\
Analysis of 18171 strong reflections:
Potential blank images: 41 -> 45
Analysis of 18171 indexed reflections:
Potential blank images: 41 -> 45
Analysis of 44340 integrated reflections:
Potential blank images: 1 -> 9
Potential blank images: 41 -> 45
"""
        in captured.out
    )
    json_file = run_in_tmp_path / "detect_blanks.json"
    assert json_file.is_file()
    with json_file.open() as fh:
        d = json.load(fh)
    assert d["strong"]["blank_regions"] == [[40, 45]]
    assert d["indexed"]["blank_regions"] == [[40, 45]]
    assert d["integrated"]["blank_regions"] == [[0, 9], [40, 45]]


def test_passing_still_images_raises_sysexit(dials_data, run_in_tmp_path):
    path = dials_data("thaumatin_grid_scan", pathlib=True)

    # Import the data
    result = procrunner.run(
        ["dials.import", "output.experiments=imported.expt"]
        + list(path.glob("*.cbf*")),
        working_directory=run_in_tmp_path,
    )
    assert not result.returncode and not result.stderr
    expts_file = run_in_tmp_path / "imported.expt"
    assert expts_file.is_file()

    # Find the spots
    result = procrunner.run(
        ["dials.find_spots", "imported.expt", "min_spot_size=3", "nproc=1"],
        working_directory=run_in_tmp_path,
    )
    assert not result.returncode and not result.stderr
    refl_file = run_in_tmp_path / "strong.refl"
    assert refl_file.is_file()

    # Now test that passing stills raises a sys.exit
    with pytest.raises(SystemExit):
        detect_blanks.run([str(expts_file), str(refl_file)])
