from __future__ import absolute_import, division, print_function

import procrunner


def test_datablock(dials_data, tmpdir):
    input_filename = dials_data("centroid_test_data").join("datablock.json").strpath
    mask_filename = dials_data("centroid_test_data").join("lookup_mask.pickle").strpath
    output_filename = tmpdir.join("output_experiments.json").strpath

    result = procrunner.run(
        [
            "dials.apply_mask",
            "input.datablock=%s" % input_filename,
            "input.mask=%s" % mask_filename,
            "output.datablock=%s" % output_filename,
        ],
        working_directory=tmpdir.strpath,
    )
    assert result["exitcode"] == 0
    assert result["stderr"] == ""

    from dials.array_family import flex  # noqa: F401, import dependency
    from dxtbx.datablock import DataBlockFactory

    datablocks = DataBlockFactory.from_json_file(output_filename)

    assert len(datablocks) == 1
    imagesets = datablocks[0].extract_imagesets()
    assert len(imagesets) == 1
    imageset = imagesets[0]
    assert imageset.external_lookup.mask.filename == mask_filename


def test_experiments(dials_data, tmpdir):
    input_filename = dials_data("centroid_test_data").join("experiments.json").strpath
    mask_filename = dials_data("centroid_test_data").join("lookup_mask.pickle").strpath
    output_filename = tmpdir.join("output_experiments.json").strpath

    result = procrunner.run(
        [
            "dials.apply_mask",
            "input.experiments=%s" % input_filename,
            "input.mask=%s" % mask_filename,
            "output.experiments=%s" % output_filename,
        ],
        working_directory=tmpdir.strpath,
    )
    assert result["exitcode"] == 0
    assert result["stderr"] == ""

    from dials.array_family import flex  # noqa: F401, import dependency
    from dxtbx.model.experiment_list import ExperimentListFactory

    experiments = ExperimentListFactory.from_json_file(output_filename)

    assert len(experiments) == 1
    imageset = experiments[0].imageset
    assert imageset.external_lookup.mask.filename == mask_filename
