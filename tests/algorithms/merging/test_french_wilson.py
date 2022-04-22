from __future__ import annotations

from dxtbx.serialize import load

from dials.algorithms.merging import french_wilson, merge
from dials.array_family import flex

MAX_SIGF_OVER_F_ACENTRIC = max(
    sigf / f for f, sigf in zip(french_wilson.ac_zf, french_wilson.ac_zf_sd)
)
MAX_SIGF_OVER_F_CENTRIC = max(
    sigf / f for f, sigf in zip(french_wilson.c_zf, french_wilson.c_zf_sd)
)


def test_french_wilson_insulin(dials_data):
    insulin = dials_data("insulin_processed", pathlib=True)
    expts = load.experiment_list(insulin / "scaled.expt", check_format=False)
    refls = flex.reflection_table.from_file(insulin / "scaled.refl")
    merged, _, _ = merge.merge(expts, refls)
    merged_intensities = merged.array()

    amplitudes = french_wilson.french_wilson(merged_intensities)
    assert amplitudes.is_xray_amplitude_array()
    assert flex.min(amplitudes.data()) >= 0
    assert flex.min(amplitudes.sigmas()) >= 0
    sigF_over_F = amplitudes.sigmas() / amplitudes.data()
    is_centric = amplitudes.centric_flags().data()
    assert flex.max(sigF_over_F.select(~is_centric)) <= MAX_SIGF_OVER_F_ACENTRIC
    assert flex.max(sigF_over_F.select(is_centric)) <= MAX_SIGF_OVER_F_CENTRIC
