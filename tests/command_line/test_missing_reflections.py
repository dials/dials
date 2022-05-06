from __future__ import annotations

from dials.command_line import missing_reflections


def test_l_cysteine_4_sweeps_scaled(dials_data, capsys):
    data = dials_data("l_cysteine_4_sweeps_scaled", pathlib=True)
    missing_reflections.run(
        args=[str(data / "scaled_30.expt"), str(data / "scaled_30.refl")]
    )
    captured = capsys.readouterr()
    assert "Completeness in resolution range: 0.754473" in captured.out
    assert "Completeness with d_max=infinity: 0.753543" in captured.out
    assert "# reflections |   % missing | Resolution range (Å)" in captured.out
    assert "260 |        16   | 1.37-0.59" in captured.out


def test_vmxi_proteinase_k_sweeps_integrated(dials_data, capsys):
    data = dials_data("vmxi_proteinase_k_sweeps", pathlib=True)
    missing_reflections.run(
        args=[
            str(data / "experiments_0.expt"),
            str(data / "reflections_0.refl"),
            str(data / "experiments_1.expt"),
            str(data / "reflections_1.refl"),
        ]
    )
    captured = capsys.readouterr()
    assert "Completeness in resolution range: 0.781833" in captured.out
    assert "Completeness with d_max=infinity: 0.7818" in captured.out
    assert "# reflections |   % missing | Resolution range (Å)" in captured.out
    assert "4899 |        20.7 | 2.36-1.80" in captured.out
    assert "190 |         0.8 | 2.36-1.80" in captured.out


def test_insulin_scaled(dials_data, capsys):
    data = dials_data("insulin_processed", pathlib=True)
    missing_reflections.run(args=[str(data / "scaled.expt"), str(data / "scaled.refl")])
    captured = capsys.readouterr()
    assert "Resolution range: 55.2195 1.45064" in captured.out
    assert "Completeness in resolution range: 0.792288" in captured.out
    assert "Completeness with d_max=infinity: 0.792288" in captured.out
    assert "# reflections |   % missing | Resolution range (Å)" in captured.out
    assert (
        "2925 |        20.6 | 1.84-1.45" in captured.out
        or "2924 |        20.6 | 1.84-1.45" in captured.out
    )
    assert "163 |         1.1 | 1.57-1.45" in captured.out


def test_insulin_scaled_d_min_d_max(dials_data, capsys):
    data = dials_data("insulin_processed", pathlib=True)
    missing_reflections.run(
        args=[
            str(data / "scaled.expt"),
            str(data / "scaled.refl"),
            "d_min=1.863199",  # inscribed circle
            "d_max=55",
            "min_component_size=10",
        ]
    )
    captured = capsys.readouterr()
    assert "Resolution range: 39.0461 1.86463" in captured.out
    assert "Completeness in resolution range: 0.996462" in captured.out
    assert "Completeness with d_max=infinity: 0.996315" in captured.out
    assert "No connected regions of missing reflections identified" in captured.out
