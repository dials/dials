"""Test for the preservation of identifiers in dials processing."""
import procrunner

from dxtbx.serialize import load

from dials.array_family import flex


def test_preservation_of_identifiers(dials_data, tmpdir):
    """Run the dials processing workflow, checking for preservation of identifiers.

    This is just a simple case. The individual programs that are expected to
    change the identifiers are tested separately, this is to check that the
    other programs maintain the identifiers through processing.
    """

    # First import - should set a unique id.
    image_files = dials_data("centroid_test_data").listdir("centroid*.cbf", sort=True)
    result = procrunner.run(
        ["dials.import", "output.experiments=imported.expt"]
        + [f.strpath for f in image_files],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("imported.expt").check(file=1)

    imported_exp_path = tmpdir.join("imported.expt").strpath
    experiments = load.experiment_list(imported_exp_path)
    import_expt_id = experiments[0].identifier
    assert import_expt_id != ""

    # Now find spots.
    result = procrunner.run(
        [
            "dials.find_spots",
            "nproc=1",
            imported_exp_path,
            "output.reflections=strong.refl",
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("strong.refl").check(file=1)

    strong_refl_path = tmpdir.join("strong.refl").strpath
    reflections = flex.reflection_table.from_file(strong_refl_path)
    assert dict(reflections.experiment_identifiers()) == {0: import_expt_id}

    # Now index
    result = procrunner.run(
        [
            "dials.index",
            strong_refl_path,
            imported_exp_path,
            "output.reflections=indexed.refl",
            "output.experiments=indexed.expt",
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("indexed.refl").check(file=1)
    assert tmpdir.join("indexed.expt").check(file=1)

    indexed_exp_path = tmpdir.join("indexed.expt").strpath
    experiments = load.experiment_list(indexed_exp_path)
    indexed_refl_path = tmpdir.join("indexed.refl").strpath
    reflections = flex.reflection_table.from_file(indexed_refl_path)

    indexed_expt_id = experiments[0].identifier
    assert indexed_expt_id != ""

    assert list(experiments.identifiers()) == [indexed_expt_id]
    assert dict(reflections.experiment_identifiers()) == {0: indexed_expt_id}

    # Now refine bravais setting
    result = procrunner.run(
        ["dials.refine_bravais_settings", indexed_refl_path, indexed_exp_path],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("bravais_setting_9.expt").check(file=1)

    bravais_exp_path = tmpdir.join("bravais_setting_9.expt").strpath
    experiments = load.experiment_list(bravais_exp_path)
    assert experiments[0].identifier == indexed_expt_id

    # Now reindex
    result = procrunner.run(
        [
            "dials.reindex",
            indexed_refl_path,
            indexed_exp_path,
            "change_of_basis_op=b,c,a",
            "output.reflections=reindexed.refl",
            "output.experiments=reindexed.expt",
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("reindexed.expt").check(file=1)
    assert tmpdir.join("reindexed.refl").check(file=1)

    reindexed_exp_path = tmpdir.join("reindexed.expt").strpath
    experiments = load.experiment_list(reindexed_exp_path)
    reindexed_refl_path = tmpdir.join("reindexed.refl").strpath
    reflections = flex.reflection_table.from_file(reindexed_refl_path)

    assert list(experiments.identifiers()) == [indexed_expt_id]
    assert dict(reflections.experiment_identifiers()) == {0: indexed_expt_id}

    # Now refine
    result = procrunner.run(
        [
            "dials.refine",
            reindexed_refl_path,
            reindexed_exp_path,
            "output.reflections=refined.refl",
            "output.experiments=refined.expt",
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("refined.expt").check(file=1)
    assert tmpdir.join("refined.refl").check(file=1)

    refined_exp_path = tmpdir.join("refined.expt").strpath
    experiments = load.experiment_list(refined_exp_path)
    refined_refl_path = tmpdir.join("refined.refl").strpath
    reflections = flex.reflection_table.from_file(refined_refl_path)

    assert list(experiments.identifiers()) == [indexed_expt_id]
    assert dict(reflections.experiment_identifiers()) == {0: indexed_expt_id}

    # Now integrate
    result = procrunner.run(
        [
            "dials.integrate",
            "nproc=1",
            refined_refl_path,
            refined_exp_path,
            "output.reflections=integrated.refl",
            "output.experiments=integrated.expt",
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("integrated.expt").check(file=1)
    assert tmpdir.join("integrated.refl").check(file=1)

    integrated_exp_path = tmpdir.join("integrated.expt").strpath
    experiments = load.experiment_list(integrated_exp_path)
    integrated_refl_path = tmpdir.join("integrated.refl").strpath
    reflections = flex.reflection_table.from_file(integrated_refl_path)

    assert list(experiments.identifiers()) == [indexed_expt_id]
    assert dict(reflections.experiment_identifiers()) == {0: indexed_expt_id}

    # Now run cosym (symmetry fails due to small amount of data)
    result = procrunner.run(
        [
            "dials.symmetry",
            integrated_refl_path,
            integrated_exp_path,
            "output.reflections=symmetrized.refl",
            "output.experiments=symmetrized.expt",
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("symmetrized.expt").check(file=1)
    assert tmpdir.join("symmetrized.refl").check(file=1)

    symmetrized_exp_path = tmpdir.join("symmetrized.expt").strpath
    experiments = load.experiment_list(symmetrized_exp_path)
    symmetrized_refl_path = tmpdir.join("symmetrized.refl").strpath
    reflections = flex.reflection_table.from_file(symmetrized_refl_path)

    assert list(experiments.identifiers()) == [indexed_expt_id]
    assert dict(reflections.experiment_identifiers()) == {0: indexed_expt_id}

    # Now scale
    result = procrunner.run(
        [
            "dials.scale",
            symmetrized_refl_path,
            symmetrized_exp_path,
            "output.reflections=scaled.refl",
            "output.experiments=scaled.expt",
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("scaled.expt").check(file=1)
    assert tmpdir.join("scaled.refl").check(file=1)

    scaled_exp_path = tmpdir.join("scaled.expt").strpath
    experiments = load.experiment_list(scaled_exp_path)
    scaled_refl_path = tmpdir.join("scaled.refl").strpath
    reflections = flex.reflection_table.from_file(scaled_refl_path)

    assert list(experiments.identifiers()) == [indexed_expt_id]
    assert dict(reflections.experiment_identifiers()) == {0: indexed_expt_id}

    # Now do two-theta refine
    result = procrunner.run(
        [
            "dials.two_theta_refine",
            scaled_refl_path,
            scaled_exp_path,
            "output.experiments=tt.expt",
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("tt.expt").check(file=1)

    tt_exp_path = tmpdir.join("tt.expt").strpath
    experiments = load.experiment_list(tt_exp_path)

    assert list(experiments.identifiers()) == [indexed_expt_id]
