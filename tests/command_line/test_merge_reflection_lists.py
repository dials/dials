from __future__ import annotations

import shutil
import subprocess

from dxtbx.util import ersatz_uuid4


def test(tmp_path):
    from dials.array_family import flex

    table = flex.reflection_table()
    table["hkl"] = flex.miller_index(360)
    table["id"] = flex.int(360)
    table["intensity.sum.value"] = flex.double(360)
    for id in set(table["id"]):
        table.experiment_identifiers()[id] = ersatz_uuid4()
    table.as_file(tmp_path / "temp1.refl")
    table.as_file(tmp_path / "temp2.refl")

    result = subprocess.run(
        [
            shutil.which("dials.merge_reflection_lists"),
            tmp_path / "temp1.refl",
            tmp_path / "temp2.refl",
            "method=update",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    table = flex.reflection_table.from_file(tmp_path / "merged.refl")
    assert len(table) == 360

    result = subprocess.run(
        [
            shutil.which("dials.merge_reflection_lists"),
            tmp_path / "temp1.refl",
            tmp_path / "temp2.refl",
            "method=extend",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    table = flex.reflection_table.from_file(tmp_path / "merged.refl")
    assert len(table) == 720
