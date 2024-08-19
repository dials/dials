from __future__ import annotations

import copy

import pytest

from dxtbx.serialize import load

from dials.algorithms.profile_model.ellipsoid.algorithm import run_ellipsoid_refinement
from dials.algorithms.profile_model.ellipsoid.indexer import reindex
from dials.array_family import flex


@pytest.mark.xdist_group(name="group1")
def test_run_ellipsoid_refinement(dials_data):
    # Download data set and the internally referenced images
    ssx = dials_data("cunir_serial_processed", pathlib=True)
    dials_data("cunir_serial", pathlib=True)

    refls = flex.reflection_table.from_file(ssx / "indexed.refl")
    expts = load.experiment_list(ssx / "indexed.expt", check_format=False)[0:1]
    refls = refls.select_on_experiment_identifiers([expts[0].identifier])
    refls = reindex(refls, expts[0])

    e_simple1 = {"spherical": 0.00027709}
    e_simple6 = {"min": 5.93276e-5, "mid": 0.00019659, "max": 0.00043976}

    initial_crystal = copy.deepcopy(expts[0].crystal)
    for model, expected in zip(
        ["simple1", "simple6"],
        [e_simple1, e_simple6],
    ):
        elist = copy.deepcopy(expts[0:1])
        out_expt, _, __ = run_ellipsoid_refinement(
            elist, refls, 0.00062, profile_model=model
        )
        for k, v in out_expt.profiles()[0].mosaicity().items():
            assert expected[k] == pytest.approx(v, abs=1e-4)
        assert (
            out_expt[0].crystal.get_unit_cell().parameters()
            != initial_crystal.get_unit_cell().parameters()
        )
        assert list(out_expt[0].crystal.get_A()) != list(initial_crystal.get_A())
        del elist[0].crystal.mosaicity

    # now try fix uc
    elist = copy.deepcopy(expts[0:1])
    out_expt, _, __ = run_ellipsoid_refinement(
        elist,
        refls,
        0.00062,
        profile_model="simple6",
        fix_orientation=True,
        fix_unit_cell=True,
    )
    for k, v in out_expt.profiles()[0].mosaicity().items():
        assert e_simple6[k] == pytest.approx(v, abs=1e-4)
    assert out_expt[0].crystal.get_unit_cell().parameters() == pytest.approx(
        initial_crystal.get_unit_cell().parameters(), abs=1e-12
    )
    assert list(out_expt[0].crystal.get_A()) == pytest.approx(
        list(initial_crystal.get_A()), abs=1e-12
    )
