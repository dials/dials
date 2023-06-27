from __future__ import annotations

from unittest import mock

from boilerplate import run


def test_boilerplate(dials_data, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    x4wide = dials_data("x4wide_processed", pathlib=True)
    with mock.patch("boilerplate.dials.util.log"):
        run(
            args=[
                str(x4wide / "AUTOMATIC_DEFAULT_scaled.expt"),
                str(x4wide / "AUTOMATIC_DEFAULT_scaled.refl"),
                "integer_parameter=42",
                "bool_parameter=True",
            ]
        )

    assert (tmp_path / "stronger.refl").is_file()
