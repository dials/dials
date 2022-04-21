from __future__ import annotations

import procrunner
import pytest

from iotbx import ccp4_map
from scitbx.array_family import flex


def test_rs_mapper(dials_data, tmp_path):
    result = procrunner.run(
        [
            "dials.rs_mapper",
            dials_data("centroid_test_data", pathlib=True)
            / "imported_experiments.json",
            'map_file="junk.ccp4"',
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "junk.ccp4").is_file()

    # load results
    m = ccp4_map.map_reader(file_name=str(tmp_path / "junk.ccp4"))
    assert len(m.data) == 7189057
    assert m.header_min == 0.0
    assert flex.min(m.data) == 0.0

    assert m.header_max == 2052.75
    assert flex.max(m.data) == 2052.75

    assert m.header_mean == pytest.approx(0.018924040719866753, abs=1e-6)
    assert flex.mean(m.data) == pytest.approx(0.01892407052218914, abs=1e-6)


def test_masked(dials_data, tmp_path):
    procrunner.run(
        [
            "dials.import",
            dials_data("thaumatin_eiger_screen", pathlib=True) / "Therm_6_1_master.h5",
        ],
        working_directory=tmp_path,
    )
    result = procrunner.run(
        ["dials.rs_mapper", "imported.expt", "map_file=junk.ccp4"],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "junk.ccp4").is_file()

    # load results
    m = ccp4_map.map_reader(file_name=str(tmp_path / "junk.ccp4"))
    assert m.header_max == pytest.approx(6330.33350)
