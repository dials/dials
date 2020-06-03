from dials.command_line import missing_reflections


def test_l_cysteine_4_sweeps_scaled(dials_data, capsys):
    missing_reflections.run(
        args=[
            (dials_data("l_cysteine_4_sweeps_scaled") / "scaled_30.expt").strpath,
            (dials_data("l_cysteine_4_sweeps_scaled") / "scaled_30.refl").strpath,
        ]
    )
    captured = capsys.readouterr()
    assert "260 reflections (16.0%): 1.37-0.59 Å" in captured.out


def test_vmxi_proteinase_k_sweeps_integrated(dials_data, capsys):
    missing_reflections.run(
        args=[
            (dials_data("vmxi_proteinase_k_sweeps") / "experiments_0.expt").strpath,
            (dials_data("vmxi_proteinase_k_sweeps") / "reflections_0.refl").strpath,
        ]
    )
    captured = capsys.readouterr()
    assert "6648 reflections (28.2%): 2.51-1.80 Å" in captured.out
    assert "307 reflections (1.3%): 103.98-1.80 Å" in captured.out


def test_insulin_scaled(dials_data, capsys):
    missing_reflections.run(
        args=[
            (dials_data("insulin_processed") / "scaled.expt").strpath,
            (dials_data("insulin_processed") / "scaled.refl").strpath,
        ]
    )
    captured = capsys.readouterr()
    assert "2925 reflections (20.6%): 1.84-1.45 Å" in captured.out
    assert "163 reflections (1.1%): 1.57-1.45 Å" in captured.out
