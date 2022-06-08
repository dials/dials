from __future__ import annotations

import pytest

from cctbx import sgtbx
from dxtbx.serialize import load

from dials.algorithms.merging import french_wilson, merge
from dials.array_family import flex

MAX_SIGF_OVER_F_ACENTRIC = max(
    sigf / f for f, sigf in zip(french_wilson.ac_zf, french_wilson.ac_zf_sd)
)
MAX_SIGF_OVER_F_CENTRIC = max(
    sigf / f for f, sigf in zip(french_wilson.c_zf, french_wilson.c_zf_sd)
)


def check_french_wilson_amplitudes(amplitudes):
    assert amplitudes.is_xray_amplitude_array()
    assert flex.min(amplitudes.data()) >= 0
    assert flex.min(amplitudes.sigmas()) >= 0
    sigF_over_F = amplitudes.sigmas() / amplitudes.data()
    is_centric = amplitudes.centric_flags().data()
    assert flex.max(sigF_over_F.select(~is_centric)) <= MAX_SIGF_OVER_F_ACENTRIC
    if is_centric.count(True):
        assert flex.max(sigF_over_F.select(is_centric)) <= MAX_SIGF_OVER_F_CENTRIC


def test_french_wilson_insulin(dials_data):
    insulin = dials_data("insulin_processed", pathlib=True)
    expts = load.experiment_list(insulin / "scaled.expt", check_format=False)
    refls = flex.reflection_table.from_file(insulin / "scaled.refl")
    merged, _, _ = merge.merge(expts, refls)
    merged_intensities = merged.array()

    amplitudes = french_wilson.french_wilson(merged_intensities)
    check_french_wilson_amplitudes(amplitudes)


@pytest.mark.parametrize("space_group_symbol", [None, "P1"])
def test_french_wilson_l_cysteine(dials_data, space_group_symbol):
    l_cysteine = dials_data("l_cysteine_4_sweeps_scaled")
    expts = load.experiment_list(l_cysteine / "scaled_30.expt", check_format=False)
    refls = flex.reflection_table.from_file(l_cysteine / "scaled_30.refl")
    if space_group_symbol:
        expts[0].crystal.set_space_group(
            sgtbx.space_group_info(space_group_symbol).group()
        )
    merged, _, _ = merge.merge(expts, refls)
    merged_intensities = merged.array()

    amplitudes = french_wilson.french_wilson(merged_intensities)
    check_french_wilson_amplitudes(amplitudes)
