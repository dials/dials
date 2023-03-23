from __future__ import annotations

import pytest

from cctbx import crystal, miller, sgtbx
from scitbx.array_family import flex

from dials.algorithms.merging.merge import (
    MTZDataClass,
    generate_r_free_flags,
    make_merged_mtz_file,
    r_free_flags_from_reference,
)
from dials.algorithms.merging.reporting import dano_over_sigdano
from dials.command_line.merge import phil_scope


def test_dano_over_sigdano():
    """Create anomalous difference data and check the calculated value."""
    ms = miller.build_set(
        crystal_symmetry=crystal.symmetry(
            space_group_symbol="P222", unit_cell=(6, 6, 6, 90, 90, 90)
        ),
        anomalous_flag=True,
        d_min=5.0,
    ).expand_to_p1()
    ma = miller.array(
        ms, data=flex.double([1, 2, 1, 3, 1, 4]), sigmas=flex.double(6, 1)
    )
    # differences are (1, 2, 3) i.e. mean 2, sigmas (sqrt2, sqrt2, sqrt2)
    assert dano_over_sigdano(ma) == pytest.approx(2**0.5)


def test_generate_r_free_flags():
    ms = miller.build_set(
        crystal_symmetry=sgtbx.space_group_info("P422").any_compatible_crystal_symmetry(
            volume=1e5
        ),
        anomalous_flag=False,
        d_min=2,
    )
    ma = miller.array(
        ms, data=flex.double(ms.size(), 1), sigmas=flex.double(ms.size(), 1)
    )
    ma = ma.select(flex.random_bool(ma.size(), 0.9))
    mtz_datasets = [MTZDataClass(merged_array=ma)]
    params = phil_scope.extract()
    r_free_flags = generate_r_free_flags(params, mtz_datasets)
    assert (
        pytest.approx(
            (r_free_flags.data() == 0).count(True) / r_free_flags.size(), abs=0.01
        )
        == 0.05
    )
    assert set(r_free_flags.data()) == set(range(20))
    assert pytest.approx(r_free_flags.completeness(), rel=0.02) == 0.9

    params.r_free_flags.relative_to_complete_set = True
    r_free_flags_complete = generate_r_free_flags(params, mtz_datasets)
    assert r_free_flags_complete.completeness() == 1
    assert r_free_flags_complete.size() > r_free_flags.size()
    assert r_free_flags_complete.d_min() == pytest.approx(2, rel=1e-4)

    params.r_free_flags.d_min = 1.5
    r_free_flags_d_min = generate_r_free_flags(params, mtz_datasets)
    assert r_free_flags_d_min.d_min() == pytest.approx(1.5, rel=1e-4)


def test_r_free_flags_from_reference(tmp_path):
    # First generate some r-free flags and save to free.mtz
    mtz_file = tmp_path / "free.mtz"
    ms = miller.build_set(
        crystal_symmetry=sgtbx.space_group_info("P422").any_compatible_crystal_symmetry(
            volume=1e5
        ),
        anomalous_flag=False,
        d_min=2,
    )
    ma = miller.array(
        ms, data=flex.double(ms.size(), 1), sigmas=flex.double(ms.size(), 1)
    )
    ma = ma.select(flex.random_bool(ma.size(), 0.9))
    mtz_datasets = [MTZDataClass(merged_array=ma)]
    params = phil_scope.extract()
    r_free_flags = generate_r_free_flags(params, mtz_datasets)
    mtz = make_merged_mtz_file(mtz_datasets, r_free_array=r_free_flags)
    mtz.write(str(mtz_file))

    # Now actually test r_free_flags_from_reference
    params.r_free_flags.reference = str(mtz_file)
    r_free_flags = r_free_flags_from_reference(params, mtz_datasets)
    assert pytest.approx(r_free_flags.completeness(), rel=0.02) == 0.9

    # And now with extension of flags to complete set
    params.r_free_flags.extend = True
    r_free_flags_extended = r_free_flags_from_reference(params, mtz_datasets)
    assert pytest.approx(r_free_flags_extended.completeness(), rel=5e-4) == 1
