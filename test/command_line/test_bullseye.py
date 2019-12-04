from __future__ import absolute_import, division, print_function

from dxtbx.model.experiment_list import ExperimentListFactory
from dials.array_family import flex


def test_bullseye_select_scans_from_reflections_1(dials_data, tmpdir):
    location = dials_data("i19_1_pdteet_index")
    refls = location.join("indexed.refl")
    expts = location.join("indexed.expt")

    data = flex.reflection_table.from_file(refls.strpath)

    # prune only indexed reflections

    data = data.select(data.get_flags(data.flags.indexed))
    z = data["xyzobs.px.value"].parts()[2]

    sel = ((data["id"] == 0) & (z > 200)) | ((data["id"] == 1) & (z < 650))
    trimmed = data.select(sel)

    from dials.command_line.bullseye import select_scans_from_reflections

    expts = ExperimentListFactory.from_json_file(expts.strpath, check_format=False)

    for j, e in enumerate(expts):
        sel = trimmed.select(trimmed["id"] == j)
        scans = select_scans_from_reflections(sel, e.scan)
        assert len(scans) == 1
        if j == 0:
            assert scans[0].get_image_range()[1] == e.scan.get_image_range()[1]
        else:
            assert scans[0].get_image_range()[0] == e.scan.get_image_range()[0]


def test_bullseye_select_scans_from_reflections_2(dials_data, tmpdir):
    location = dials_data("i19_1_pdteet_index")
    refls = location.join("indexed.refl")
    expts = location.join("indexed.expt")

    data = flex.reflection_table.from_file(refls.strpath)

    data = data.select(data.get_flags(data.flags.indexed))
    z = data["xyzobs.px.value"].parts()[2]

    sel = ((data["id"] == 0) & (z < 200)) | ((data["id"] == 0) & (z > 650))
    trimmed = data.select(sel)

    from dials.command_line.bullseye import select_scans_from_reflections

    e = ExperimentListFactory.from_json_file(expts.strpath, check_format=False)[0]

    scans = select_scans_from_reflections(trimmed, e.scan)
    assert len(scans) == 2
    assert scans[0].get_image_range()[0] == e.scan.get_image_range()[0]
    assert scans[1].get_image_range()[1] == e.scan.get_image_range()[1]
