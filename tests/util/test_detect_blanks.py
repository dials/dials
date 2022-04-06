from __future__ import annotations

import pytest

from dxtbx.model import ExperimentList

from dials.array_family import flex
from dials.util import detect_blanks


def test_blank_counts_analysis(dials_data):
    expts = ExperimentList.from_file(
        dials_data("insulin_processed", pathlib=True) / "imported.expt",
        check_format=False,
    )
    refl = flex.reflection_table.from_file(
        dials_data("insulin_processed", pathlib=True) / "strong.refl"
    )
    results = detect_blanks.blank_counts_analysis(
        refl, expts[0].scan, phi_step=5, fractional_loss=0.1
    )
    assert set(results) == {"data", "layout", "blank_regions"}
    assert results["data"][0]["x"] == [
        2.5,
        7.5,
        12.5,
        17.5,
        22.5,
        27.5,
        32.5,
        37.5,
        42.5,
    ]
    assert results["data"][0]["y"] == [
        2827,
        2589,
        2502,
        2464,
        2515,
        2490,
        2441,
        2500,
        2505,
    ]
    assert not any(results["data"][0]["blank"])
    assert results["blank_regions"] == []

    z = refl["xyzobs.px.value"].parts()[2]
    refl_subset = refl.select((z < 10) | (z > 20))
    results = detect_blanks.blank_counts_analysis(
        refl_subset, expts[0].scan, phi_step=2, fractional_loss=0.1
    )
    assert results["data"][0]["blank"].count(True) == 5
    assert results["blank_regions"] == [(10, 21)]


def test_blank_integrated_analysis(dials_data):
    expts = ExperimentList.from_file(
        dials_data("insulin_processed", pathlib=True) / "integrated.expt",
        check_format=False,
    )
    refl = flex.reflection_table.from_file(
        dials_data("insulin_processed", pathlib=True) / "integrated.refl"
    )
    results = detect_blanks.blank_integrated_analysis(
        refl, expts[0].scan, phi_step=5, fractional_loss=0.1
    )
    assert results["data"][0]["x"] == [
        2.5,
        7.5,
        12.5,
        17.5,
        22.5,
        27.5,
        32.5,
        37.5,
        42.5,
    ]
    assert results["data"][0]["y"] == pytest.approx(
        [
            27.903266149430973,
            25.832527090455052,
            26.9236206883069,
            26.50234804728626,
            26.41019377727383,
            25.810676090828185,
            24.844906790823064,
            25.89992001081651,
            25.580718362291474,
        ]
    )
    assert not any(results["data"][0]["blank"])
    assert results["blank_regions"] == []

    # Now with some "blank" regions - make some of the reflections weak
    z = refl["xyzobs.px.value"].parts()[2]
    refl["intensity.prf.value"].set_selected(z < 10, refl["intensity.prf.value"] * 0.05)
    results = detect_blanks.blank_integrated_analysis(
        refl, expts[0].scan, phi_step=5, fractional_loss=0.1
    )
    assert results["data"][0]["y"] == pytest.approx(
        [
            1.3951633074715482,
            1.2916263545227527,
            26.9236206883069,
            26.50234804728626,
            26.41019377727383,
            25.810676090828185,
            24.844906790823064,
            25.89992001081651,
            25.580718362291474,
        ]
    )
    assert results["data"][0]["blank"] == [
        True,
        True,
        False,
        False,
        False,
        False,
        False,
        False,
        False,
    ]
    assert results["blank_regions"] == [(0, 10)]

    # Unset the integrated_prf flags, so the analysis should instead use the umodified
    # intensity.sum.value instead
    refl.unset_flags(flex.bool(len(refl), True), refl.flags.integrated_prf)
    results = detect_blanks.blank_integrated_analysis(
        refl, expts[0].scan, phi_step=5, fractional_loss=0.1
    )
    assert not any(results["data"][0]["blank"])
    assert results["blank_regions"] == []
