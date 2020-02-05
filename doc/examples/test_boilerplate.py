from boilerplate import run


def test_boilerplate(dials_data, tmpdir):
    x4wide_dir = dials_data("x4wide_processed")

    with tmpdir.as_cwd():
        run(
            args=[
                x4wide_dir.join("AUTOMATIC_DEFAULT_scaled.expt").strpath,
                x4wide_dir.join("AUTOMATIC_DEFAULT_scaled.refl").strpath,
                "integer_parameter=42",
                "bool_parameter=True",
            ]
        )

    assert tmpdir.join("stronger.refl").check()
