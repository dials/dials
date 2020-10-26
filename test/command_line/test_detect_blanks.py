import json

import procrunner
import pytest

from dials.array_family import flex
from dials.command_line import detect_blanks


def test_strong(dials_data, capsys, run_in_tmpdir):
    expts_file = dials_data("insulin_processed") / "imported.expt"
    refl_file = dials_data("insulin_processed") / "strong.refl"
    refl = flex.reflection_table.from_file(refl_file)
    z = refl["xyzobs.px.value"].parts()[2]
    refl_subset = refl.select((z < 10) | (z > 20))
    refl_subset.as_file("mod.refl")
    detect_blanks.run([expts_file.strpath, "mod.refl"])
    captured = capsys.readouterr()
    assert "Potential blank images: 11 -> 21" in captured.out
    json_file = run_in_tmpdir / "blanks.json"
    assert json_file.check(file=1)
    with json_file.open() as fh:
        d = json.load(fh)
        assert d["strong"]["blank_regions"] == [[10, 21]]
        assert not d["indexed"]
        assert not d["integrated"]


def test_integrated(dials_data, capsys, run_in_tmpdir):
    expts_file = dials_data("insulin_processed") / "integrated.expt"
    refl_file = dials_data("insulin_processed") / "integrated.refl"
    refl = flex.reflection_table.from_file(refl_file)
    z = refl["xyzobs.px.value"].parts()[2]
    sel = z < 40
    refl_subset = refl.select(sel)
    refl_subset["intensity.prf.value"].set_selected(
        z.select(sel) < 10, refl_subset["intensity.prf.value"] * 0.05
    )
    refl_subset.as_file("mod.refl")
    detect_blanks.run(
        [expts_file.strpath, "mod.refl", "output.json=detect_blanks.json"]
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
    json_file = run_in_tmpdir / "detect_blanks.json"
    assert json_file.check(file=1)
    with json_file.open() as fh:
        d = json.load(fh)
        assert d["strong"]["blank_regions"] == [[40, 45]]
        assert d["indexed"]["blank_regions"] == [[40, 45]]
        assert d["integrated"]["blank_regions"] == [[0, 9], [40, 45]]


def test_still_raises_sysexit(dials_data, run_in_tmpdir):
    path = dials_data("thaumatin_grid_scan")

    # Import the data
    result = procrunner.run(
        ["dials.import", "output.experiments=imported.expt"] + path.listdir("*.cbf*")
    )
    assert not result.returncode and not result.stderr
    expts_file = run_in_tmpdir.join("imported.expt")
    assert expts_file.check(file=1)

    # Find the spots
    result = procrunner.run(
        ["dials.find_spots", "imported.expt", "min_spot_size=3"],
    )
    assert not result.returncode and not result.stderr
    refl_file = run_in_tmpdir.join("strong.refl")
    assert refl_file.check(file=1)

    # Now test that passing stills raises a sys.exit
    with pytest.raises(SystemExit):
        detect_blanks.run([expts_file.strpath, refl_file.strpath])
