from __future__ import absolute_import, division, print_function

import procrunner


def test(dials_data, tmpdir):
    input_filename = dials_data("centroid_test_data").join("datablock.json").strpath
    mask_filename = dials_data("centroid_test_data").join("lookup_mask.pickle").strpath
    output_filename = tmpdir.join("output.expt").strpath

    result = procrunner.run(
        [
            "dials.apply_mask",
            "input.experiments=%s" % input_filename,
            "input.mask=%s" % mask_filename,
            "output.experiments=%s" % output_filename,
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr

    from dials.array_family import flex  # noqa: F401, import dependency
    from dxtbx.model.experiment_list import ExperimentListFactory

    experiments = ExperimentListFactory.from_json_file(output_filename)

    assert len(experiments) == 1
    imagesets = experiments.imagesets()
    assert len(imagesets) == 1
    imageset = imagesets[0]
    assert imageset.external_lookup.mask.filename == mask_filename


def test_experiments(dials_data, tmpdir):
    input_filename = dials_data("centroid_test_data").join("experiments.json").strpath
    mask_filename = dials_data("centroid_test_data").join("lookup_mask.pickle").strpath
    output_filename = tmpdir.join("output.expt").strpath

    result = procrunner.run(
        [
            "dials.apply_mask",
            "input.experiments=%s" % input_filename,
            "input.mask=%s" % mask_filename,
            "output.experiments=%s" % output_filename,
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr

    from dials.array_family import flex  # noqa: F401, import dependency
    from dxtbx.model.experiment_list import ExperimentListFactory

    experiments = ExperimentListFactory.from_json_file(output_filename)

    assert len(experiments) == 1
    imageset = experiments[0].imageset
    assert imageset.external_lookup.mask.filename == mask_filename
