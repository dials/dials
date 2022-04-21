from __future__ import annotations

import procrunner

from dxtbx.serialize import load

from dials.array_family import flex


def test_for_preservation_of_identifiers_in_dials_processing(dials_data, tmp_path):
    """Run the dials processing workflow, checking for preservation of identifiers.

    This is just a simple case. The individual programs that are expected to
    change the identifiers are tested separately, this is to check that the
    other programs maintain the identifiers through processing.
    """

    imported = tmp_path / "imported.expt"
    strong_refl = tmp_path / "strong.refl"
    indexed_expt = tmp_path / "indexed.expt"
    indexed_refl = tmp_path / "indexed.refl"
    bravais_expt = tmp_path / "bravais_setting_9.expt"
    reindexed_expt = tmp_path / "reindexed.expt"
    reindexed_refl = tmp_path / "reindexed.refl"
    refined_expt = tmp_path / "refined.expt"
    refined_refl = tmp_path / "refined.refl"
    integrated_expt = tmp_path / "integrated.expt"
    integrated_refl = tmp_path / "integrated.refl"
    symmetrized_expt = tmp_path / "symmetrized.expt"
    symmetrized_refl = tmp_path / "symmetrized.refl"
    scaled_expt = tmp_path / "scaled.expt"
    scaled_refl = tmp_path / "scaled.refl"
    tt_expt = tmp_path / "tt.expt"

    # First import - should set a unique id.
    image_files = dials_data("centroid_test_data", pathlib=True).glob("centroid*.cbf")
    result = procrunner.run(
        ["dials.import", f"output.experiments={imported.name}"] + sorted(image_files),
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr

    assert imported.is_file()
    experiments = load.experiment_list(imported)
    import_expt_id = experiments[0].identifier
    assert import_expt_id != ""

    # Now find spots.
    result = procrunner.run(
        [
            "dials.find_spots",
            "nproc=1",
            imported,
            f"output.reflections={strong_refl.name}",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr

    assert strong_refl.is_file()
    reflections = flex.reflection_table.from_file(strong_refl)
    assert dict(reflections.experiment_identifiers()) == {0: import_expt_id}

    # Now index
    result = procrunner.run(
        [
            "dials.index",
            strong_refl,
            imported,
            f"output.reflections={indexed_refl.name}",
            f"output.experiments={indexed_expt.name}",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr

    assert indexed_expt.is_file()
    assert indexed_refl.is_file()
    experiments = load.experiment_list(indexed_expt)
    reflections = flex.reflection_table.from_file(indexed_refl)

    indexed_expt_id = experiments[0].identifier
    assert indexed_expt_id
    assert list(experiments.identifiers()) == [indexed_expt_id]
    assert dict(reflections.experiment_identifiers()) == {0: indexed_expt_id}

    # Now refine bravais setting
    result = procrunner.run(
        ["dials.refine_bravais_settings", indexed_refl, indexed_expt],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr

    assert bravais_expt.is_file()
    experiments = load.experiment_list(bravais_expt)
    assert experiments[0].identifier == indexed_expt_id

    # Now reindex
    result = procrunner.run(
        [
            "dials.reindex",
            indexed_refl,
            indexed_expt,
            "change_of_basis_op=b,c,a",
            f"output.reflections={reindexed_refl.name}",
            f"output.experiments={reindexed_expt.name}",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr

    assert reindexed_expt.is_file()
    assert reindexed_refl.is_file()
    experiments = load.experiment_list(reindexed_expt)
    reflections = flex.reflection_table.from_file(reindexed_refl)

    assert list(experiments.identifiers()) == [indexed_expt_id]
    assert dict(reflections.experiment_identifiers()) == {0: indexed_expt_id}

    # Now refine
    result = procrunner.run(
        [
            "dials.refine",
            reindexed_refl,
            reindexed_expt,
            f"output.reflections={refined_refl.name}",
            f"output.experiments={refined_expt.name}",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr

    assert refined_expt.is_file()
    assert refined_refl.is_file()
    experiments = load.experiment_list(refined_expt)
    reflections = flex.reflection_table.from_file(refined_refl)

    assert list(experiments.identifiers()) == [indexed_expt_id]
    assert dict(reflections.experiment_identifiers()) == {0: indexed_expt_id}

    # Now integrate
    result = procrunner.run(
        [
            "dials.integrate",
            "nproc=1",
            refined_refl,
            refined_expt,
            f"output.reflections={integrated_refl.name}",
            f"output.experiments={integrated_expt.name}",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr

    assert integrated_expt.is_file()
    assert integrated_refl.is_file()
    experiments = load.experiment_list(integrated_expt)
    reflections = flex.reflection_table.from_file(integrated_refl)

    assert list(experiments.identifiers()) == [indexed_expt_id]
    assert dict(reflections.experiment_identifiers()) == {0: indexed_expt_id}

    # Now run cosym (symmetry fails due to small amount of data)
    result = procrunner.run(
        [
            "dials.symmetry",
            integrated_refl,
            integrated_expt,
            f"output.reflections={symmetrized_refl.name}",
            f"output.experiments={symmetrized_expt.name}",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr

    assert symmetrized_expt.is_file()
    assert symmetrized_refl.is_file()
    experiments = load.experiment_list(symmetrized_expt)
    reflections = flex.reflection_table.from_file(symmetrized_refl)

    assert list(experiments.identifiers()) == [indexed_expt_id]
    assert dict(reflections.experiment_identifiers()) == {0: indexed_expt_id}

    # Now scale
    result = procrunner.run(
        [
            "dials.scale",
            symmetrized_refl,
            symmetrized_expt,
            f"output.reflections={scaled_refl.name}",
            f"output.experiments={scaled_expt.name}",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr

    assert scaled_expt.is_file()
    assert scaled_refl.is_file()
    experiments = load.experiment_list(scaled_expt)
    reflections = flex.reflection_table.from_file(scaled_refl)

    assert list(experiments.identifiers()) == [indexed_expt_id]
    assert dict(reflections.experiment_identifiers()) == {0: indexed_expt_id}

    # Now do two-theta refine
    result = procrunner.run(
        [
            "dials.two_theta_refine",
            scaled_refl,
            scaled_expt,
            f"output.experiments={tt_expt.name}",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr

    assert tt_expt.is_file()
    experiments = load.experiment_list(tt_expt)

    assert list(experiments.identifiers()) == [indexed_expt_id]
