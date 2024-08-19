from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path
from unittest import mock

import h5py
import numpy as np
import pytest

from dxtbx.serialize import load

import dials.util.image_grouping
from dials.array_family import flex
from dials.util.image_grouping import (
    BlockInImageFile,
    ConstantMetadataForFile,
    ExtractedValues,
    FilePair,
    GroupingImageTemplates,
    ImageFile,
    ParsedYAML,
    RepeatInImageFile,
    _determine_groupings,
    example_yaml,
    get_grouping_handler,
    simple_template_example,
)


def test_yml_parsing(tmp_path):
    with open(tmp_path / "example.yaml", "w") as f:
        f.write(example_yaml)

    class MetaDataOverWrite(dials.util.image_grouping.MetadataInFile):
        def __init__(self, file, item):
            self.file = file
            self.item = item
            self.extracted_data = ExtractedValues(np.array([1.0]), False, False, True)

    with mock.patch(
        "dials.util.image_grouping.MetadataInFile", side_effect=MetaDataOverWrite
    ):
        parsed = ParsedYAML(tmp_path / "example.yaml")
        # check the images were found
        assert list(parsed._images.keys()) == [
            "/path/to/example_master.h5",
            "/path/to/example_2_master.h5",
        ]

        I1 = ImageFile("/path/to/example_master.h5", True, False)
        I2 = ImageFile("/path/to/example_2_master.h5", True, False)

        # check the metadata was found
        assert parsed.metadata_items["timepoint"] == {
            I1: dials.util.image_grouping.MetadataInFile(
                "/path/to/example_master.h5", "/timepoint"
            ),
            I2: dials.util.image_grouping.MetadataInFile(
                "/path/to/meta.h5", "/timepoint"
            ),
        }
        assert parsed.metadata_items["wavelength"] == {
            I1: ConstantMetadataForFile(0.4),
            I2: ConstantMetadataForFile(0.6),
        }

        # Now check the two groupings
        assert list(parsed.groupings.keys()) == ["merge_by", "index_by"]
        merge_by = parsed.groupings["merge_by"]
        assert merge_by.metadata_names == {"timepoint", "wavelength"}
        assert merge_by._images_to_metadata[I1] == {
            "timepoint": dials.util.image_grouping.MetadataInFile(
                "/path/to/example_master.h5", "/timepoint"
            ),
            "wavelength": ConstantMetadataForFile(0.4),
        }
        assert merge_by._images_to_metadata[I2] == {
            "timepoint": dials.util.image_grouping.MetadataInFile(
                "/path/to/meta.h5", "/timepoint"
            ),
            "wavelength": ConstantMetadataForFile(0.6),
        }
        assert merge_by.tolerances == {"timepoint": 0.1, "wavelength": 0.01}
        print(merge_by)  # test the __str__ method
        index_by = parsed.groupings["index_by"]
        assert index_by.metadata_names == {"wavelength"}
        assert index_by._images_to_metadata[I1] == {
            "wavelength": ConstantMetadataForFile(0.4)
        }
        assert index_by._images_to_metadata[I2] == {
            "wavelength": ConstantMetadataForFile(0.6)
        }
        assert index_by.tolerances == {"wavelength": 0.01}
        print(index_by)  # test the __str__ method


def test_determine_groupings(tmp_path):
    """Test the grouping with tolerances on a more involved example"""
    test_h5 = str(os.fspath(tmp_path / "meta.h5"))
    wl_example = f"""
metadata:
  timepoint:
    '/path/to/example_master.h5' : '{test_h5}:/timepoint'
  wavelength:
    '/path/to/example_master.h5' : '{test_h5}:/wavelength'
grouping:
  merge_by:
    values:
      - wavelength
      - timepoint
    tolerances:
      - 0.1
      - 0
"""
    with open(tmp_path / "example.yaml", "w") as f:
        f.write(wl_example)

    # Write some example data into a h5 file
    wl_array = np.array([1.0, 1.02, 1.05, 2.0, 2.01, 2.02])
    tp_array = np.array([1, 0, 1, 1, 1, 1], dtype=int)
    f = h5py.File(test_h5, "w")
    f.create_dataset("wavelength", data=wl_array)
    f.create_dataset("timepoint", data=tp_array)
    f.close()

    # expect three groups
    # - tp=0, wl=1.0-1.1
    # - tp=1 wl=1.0-1.1
    # - tp=1, wl=2.0-2.1

    parsed = ParsedYAML(tmp_path / "example.yaml")
    groups = _determine_groupings(parsed.groupings["merge_by"])
    assert len(groups) == 3
    # NB important that there is not an empty fourth group with tp=0, wl=1.0-1.1

    assert groups[0].min_max_for_metadata("timepoint") == (0.0, 0.0)
    assert groups[0].min_max_for_metadata("wavelength") == (1.0, 1.1)
    assert groups[1].min_max_for_metadata("timepoint") == (1.0, 1.0)
    assert groups[1].min_max_for_metadata("wavelength") == (1.0, 1.1)
    assert groups[2].min_max_for_metadata("timepoint") == (1.0, 1.0)
    assert groups[2].min_max_for_metadata("wavelength") == (2.0, 2.1)


def test_yml_parsing_template(tmp_path):
    with open(tmp_path / "example.yaml", "w") as f:
        f.write(simple_template_example)
    parsed = ParsedYAML(tmp_path / "example.yaml")

    I1 = ImageFile("/path/to/example_#####.cbf", False, True)
    # check the images were found
    assert list(parsed._images.keys()) == [I1.name]

    # check the metadata was found
    assert parsed.metadata_items["dose_point"] == {I1: RepeatInImageFile(10)}
    assert parsed.metadata_items["wavelength"] == {I1: ConstantMetadataForFile(1)}

    # Now check the two groupings
    assert list(parsed.groupings.keys()) == ["merge_by"]
    merge_by = parsed.groupings["merge_by"]
    assert merge_by.metadata_names == {"dose_point", "wavelength"}
    assert merge_by._images_to_metadata[I1] == {
        "dose_point": RepeatInImageFile(10),
        "wavelength": ConstantMetadataForFile(1),
    }
    assert merge_by.tolerances == {"dose_point": 0.1, "wavelength": 0.01}
    print(merge_by)  # test the __str__ method

    groups = _determine_groupings(merge_by)
    ftg = GroupingImageTemplates._files_to_groups(merge_by.extract_data(), groups)
    iitgi = ftg[I1]["img_idx_to_group_id"]
    for i in range(100):
        assert iitgi[i] == i % 10
    assert iitgi.single_return_val is None
    assert iitgi.group_ids.size() == 0

    simple_block_example = """
metadata:
  crystal_id:
    '/path/to/example_#####.cbf' : 'block=1:100:10'
grouping:
  merge_by:
    values:
      - crystal_id
"""
    with open(tmp_path / "block_example.yaml", "w") as f:
        f.write(simple_block_example)
    parsed = ParsedYAML(tmp_path / "block_example.yaml")

    I1 = ImageFile("/path/to/example_#####.cbf", False, True)
    # check the images were found
    assert list(parsed._images.keys()) == [I1.name]
    # check the metadata was found
    assert parsed.metadata_items["crystal_id"] == {I1: BlockInImageFile(1, 100, 10)}


invalid_example = """
---
metadata:
  timepoint:
    '/path/to/example_master.h5' : '/path/to/meta.h5:/timepoint'
  wavelength:
    '/path/to/example_master.h5' : 'repeat=2'
grouping:
  merge_by:
    values:
      - timepoint
      - wavelength
    tolerances:
      - 0.1
      - 0.01
"""
invalid_example_2 = """
---
metadata:
  timepoint:
    '/path/to/example_master.h5' : 'repeat=4'
  wavelength:
    '/path/to/example_master.h5' : 'repeat=2'
grouping:
  merge_by:
    values:
      - timepoint
      - wavelength
    tolerances:
      - 0.1
      - 0.01
"""


def test_invalid_yml(tmp_path):
    with open(tmp_path / "example.yaml", "w") as f:
        f.write(invalid_example)
    with pytest.raises(ValueError):
        _ = ParsedYAML(tmp_path / "example.yaml")
    with open(tmp_path / "example_2.yaml", "w") as f:
        f.write(invalid_example_2)
    with pytest.raises(AssertionError):
        _ = ParsedYAML(tmp_path / "example_2.yaml")


@pytest.mark.xfail(
    os.name == "nt",
    reason="Failures due to translated paths; see https://github.com/cctbx/dxtbx/issues/613",
)
def test_real_h5_example(tmp_path, dials_data):

    """This test tests a few use cases on processed data derived from h5 format."""
    fpath1 = (
        "/dls/mx/data/nt30330/nt30330-15/VMXi-AB1698/well_42/images/image_58766.nxs"
    )
    fpath2 = (
        "/dls/mx/data/nt30330/nt30330-15/VMXi-AB1698/well_39/images/image_58763.nxs"
    )
    real_example = f"""
---
metadata:
  timepoint:
    {fpath1} : 'repeat=2'
    {fpath2} : 'repeat=2'
grouping:
  group_by:
    values:
      - timepoint
"""
    # single file indices for the first dataset are 5051,5056,5058,5062,5063,5064,5065,5066
    # 5073,5074,5141,5142,5143,5144,5151,5152,5231,5248,5309
    ids_group1_file1 = [1, 2, 3, 5, 7, 9, 11, 13, 15, 17]
    ids_group2_file1 = [0, 4, 6, 8, 10, 12, 14, 16, 18]
    expected_group1_file1 = [5056, 5058, 5062, 5064, 5066, 5074, 5142, 5144, 5152, 5248]
    expected_group2_file1 = [5051, 5063, 5065, 5073, 5141, 5143, 5151, 5231, 5309]
    expected_group1_file2 = [11000, 11100, 11256, 11258, 11360, 11384, 11598]
    expected_group2_file2 = [
        11083,
        11101,
        11257,
        11361,
        11383,
        11385,
        11515,
        11599,
        11799,
    ]
    with open(tmp_path / "real_example.yaml", "w") as f:
        f.write(real_example)

    parsed = ParsedYAML(tmp_path / "real_example.yaml")
    handler = get_grouping_handler(parsed, "group_by")
    dtbp = dials_data("dtpb_serial_processed", pathlib=True)

    fps = [
        FilePair(
            dtbp / "well42_batch6_integrated.expt",
            dtbp / "well42_batch6_integrated.refl",
        )
    ]
    fd = handler.split_files_to_groups(tmp_path, fps)
    assert set(fd.keys()) == {"group_1", "group_2"}
    expts1 = load.experiment_list(fd["group_1"][0].expt, check_format=False)
    indices1 = [expt.imageset.indices()[0] for expt in expts1]
    assert indices1 == expected_group1_file1
    expts2 = load.experiment_list(fd["group_2"][0].expt, check_format=False)
    indices2 = [expt.imageset.indices()[0] for expt in expts2]
    assert indices2 == expected_group2_file1

    # Check writing the group ids to the file. Don't overwrite dials_data files though
    fps_copy = [
        FilePair(
            tmp_path / "tmp.expt",
            tmp_path / "tmp.refl",
        )
    ]
    shutil.copy(fps[0].refl, fps_copy[0].refl)
    shutil.copy(fps[0].expt, fps_copy[0].expt)
    handler.write_groupids_into_files(fps_copy)
    refls = flex.reflection_table.from_file(fps_copy[0].refl)
    assert set(refls["id"]) == set(range(19))
    sel = flex.bool(refls.size(), False)
    for id_ in ids_group1_file1:
        sel |= refls["id"] == id_
    assert set(refls["group_id"].select(sel)) == {0}
    sel = flex.bool(refls.size(), False)
    for id_ in ids_group2_file1:
        sel |= refls["id"] == id_
    assert set(refls["group_id"].select(sel)) == {1}

    # now join files to test expt files with multiple h5 images referenced:
    expts1 = load.experiment_list(
        dtbp / "well42_batch6_integrated.expt", check_format=False
    )
    expts2 = load.experiment_list(
        dtbp / "well39_batch12_integrated.expt", check_format=False
    )
    expts1.extend(expts2)
    expts1.as_file(tmp_path / "joint.expt")
    refls1 = flex.reflection_table.from_file(dtbp / "well42_batch6_integrated.refl")
    refls2 = flex.reflection_table.from_file(dtbp / "well39_batch12_integrated.refl")
    joint_refls = flex.reflection_table.concat([refls1, refls2])
    joint_refls.as_file(tmp_path / "joint.refl")

    fps = [FilePair(tmp_path / "joint.expt", tmp_path / "joint.refl")]
    fd = handler.split_files_to_groups(tmp_path, fps)
    assert set(fd.keys()) == {"group_1", "group_2"}
    expts1 = load.experiment_list(fd["group_1"][0].expt, check_format=False)
    indices1 = [expt.imageset.indices()[0] for expt in expts1]
    assert indices1 == expected_group1_file1 + expected_group1_file2
    expts2 = load.experiment_list(fd["group_2"][0].expt, check_format=False)
    indices2 = [expt.imageset.indices()[0] for expt in expts2]
    assert indices2 == expected_group2_file1 + expected_group2_file2

    test_h5 = str(os.fspath(tmp_path / "meta.h5"))

    # Write the same groupings into a h5 file
    tp_array = np.zeros((6000,), dtype=int)
    for i in expected_group1_file1:
        tp_array[i] = 0
    for i in expected_group2_file1:
        tp_array[i] = 1
    f = h5py.File(test_h5, "w")
    f.create_dataset("timepoint", data=tp_array)
    f.close()

    real_example_metafile = f"""
---
metadata:
  timepoint:
    {fpath1} : '{test_h5}:/timepoint'
    {fpath2} : 0
grouping:
  group_by:
    values:
      - timepoint
"""
    with open(tmp_path / "real_example_metafile.yaml", "w") as f:
        f.write(real_example_metafile)
    parsed = ParsedYAML(tmp_path / "real_example_metafile.yaml")
    handler = get_grouping_handler(parsed, "group_by")
    fps = [FilePair(tmp_path / "joint.expt", tmp_path / "joint.refl")]
    fd = handler.split_files_to_groups(tmp_path, fps)
    assert set(fd.keys()) == {"group_1", "group_2"}
    expts1 = load.experiment_list(fd["group_1"][0].expt, check_format=False)
    indices1 = [expt.imageset.indices()[0] for expt in expts1]
    assert indices1 == sorted(
        expected_group1_file1 + expected_group1_file2 + expected_group2_file2
    )
    expts2 = load.experiment_list(fd["group_2"][0].expt, check_format=False)
    indices2 = [expt.imageset.indices()[0] for expt in expts2]
    assert indices2 == expected_group2_file1


@pytest.mark.xfail(
    os.name == "nt",
    reason="Failures due to translated paths; see https://github.com/cctbx/dxtbx/issues/613",
)
def test_real_cbf_example(tmp_path, dials_data):

    """This test tests a few use cases on real cbf data, using the template
    metadata definition.

    First, post-indexed and pre-indexed data are split (to test the different
    imageset layouts for these kinds of files), with the "repeat=" definition.
    Then, the "block=" definition is tested on the indexed data.
    Finally, the single value definition is tested.
    """

    ssx = dials_data("cunir_serial", pathlib=True)
    fpath = str(os.fspath(ssx / "merlin0047_#####.cbf"))
    real_example = f"""
---
metadata:
  timepoint:
    '{fpath}' : 'repeat=2'
grouping:
  group_by:
    values:
      - timepoint
    tolerances:
      - 0.1
"""
    with open(tmp_path / "real_example.yaml", "w") as f:
        f.write(real_example)

    parsed = ParsedYAML(tmp_path / "real_example.yaml")
    handler = get_grouping_handler(parsed, "group_by")

    args = [shutil.which("dials.import"), f"template={fpath}"]
    result = subprocess.run(args, cwd=tmp_path, capture_output=True)
    assert not result.returncode and not result.stderr

    args = [shutil.which("dials.find_spots"), "imported.expt"]
    result = subprocess.run(args, cwd=tmp_path, capture_output=True)
    assert not result.returncode and not result.stderr
    args = [
        shutil.which("dials.ssx_index"),
        "imported.expt",
        "strong.refl",
        "unit_cell=96.4,96.4,96.4,90,90,90",
        "space_group=P213",
        "nproc=1",
        "max_lattices=2",
    ]
    result = subprocess.run(args, cwd=tmp_path, capture_output=True)
    assert not result.returncode and not result.stderr

    # First test on indexed data - here each image has its own imageset
    # only images 17001, 17002, 17003, 17004 get indexed, so expect these to be split into groups [1,0,1,0]
    fps = [FilePair(Path(tmp_path / "indexed.expt"), Path(tmp_path / "indexed.refl"))]
    fd = handler.split_files_to_groups(tmp_path, fps)

    # with max lattices=2, 17001 has two lattices, 17002,17003,17004 have one
    assert list(fd.keys()) == ["group_1", "group_2"]
    filelist_1 = fd["group_1"]
    assert len(filelist_1) == 1
    expts1 = load.experiment_list(filelist_1[0].expt)
    assert len(expts1) == 2
    assert expts1[0].scan.get_image_range() == (17002, 17002)
    assert expts1[1].scan.get_image_range() == (17004, 17004)
    filelist_2 = fd["group_2"]
    assert len(filelist_2) == 1
    expts2 = load.experiment_list(filelist_2[0].expt)
    assert len(expts2) == 3
    assert expts2[0].scan.get_image_range() == (17001, 17001)
    assert expts2[1].scan.get_image_range() == (17001, 17001)
    assert expts2[2].scan.get_image_range() == (17003, 17003)

    # Now test on imported data. Here, we have one imagesequence, expect
    # images 17000-17004, to be split into alternating groups.
    fps = [FilePair(Path(tmp_path / "imported.expt"), Path(tmp_path / "strong.refl"))]
    fd = handler.split_files_to_groups(tmp_path, fps)
    assert list(fd.keys()) == ["group_1", "group_2"]
    filelist_1 = fd["group_1"]
    assert len(filelist_1) == 1
    expts1 = load.experiment_list(filelist_1[0].expt)
    assert len(expts1) == 3
    assert expts1[0].scan.get_image_range()[0] == 17000
    assert expts1[1].scan.get_image_range()[0] == 17002
    assert expts1[2].scan.get_image_range()[0] == 17004
    filelist_2 = fd["group_2"]
    assert len(filelist_2) == 1
    expts2 = load.experiment_list(filelist_2[0].expt)
    assert len(expts2) == 2
    assert expts2[0].scan.get_image_range()[0] == 17001
    assert expts2[1].scan.get_image_range()[0] == 17003

    real_example_block = f"""
---
metadata:
  timepoint:
    '{fpath}' : 'block=17000:17004:2'
grouping:
  group_by:
    values:
      - timepoint
    tolerances:
      - 0.1
"""

    with open(tmp_path / "real_example_block.yaml", "w") as f:
        f.write(real_example_block)

    parsed = ParsedYAML(tmp_path / "real_example_block.yaml")
    handler = get_grouping_handler(parsed, "group_by")

    fps = [FilePair(Path(tmp_path / "indexed.expt"), Path(tmp_path / "indexed.refl"))]
    fd = handler.split_files_to_groups(tmp_path, fps)
    assert list(fd.keys()) == ["group_1", "group_2", "group_3"]
    # with max lattices=2, 17001 has two lattices, 17002,17003,17004 have one
    filelist_1 = fd["group_1"]
    assert len(filelist_1) == 1
    expts1 = load.experiment_list(filelist_1[0].expt)
    assert len(expts1) == 2
    assert expts1[0].scan.get_image_range()[0] == 17001
    assert expts1[1].scan.get_image_range()[0] == 17001
    filelist_2 = fd["group_2"]
    assert len(filelist_2) == 1
    expts2 = load.experiment_list(filelist_2[0].expt)
    assert len(expts2) == 2
    assert expts2[0].scan.get_image_range()[0] == 17002
    assert expts2[1].scan.get_image_range()[0] == 17003
    filelist_3 = fd["group_3"]
    assert len(filelist_3) == 1
    expts3 = load.experiment_list(filelist_3[0].expt)
    assert len(expts3) == 1
    assert expts3[0].scan.get_image_range()[0] == 17004

    # Check writing the group ids to the file. Don't overwrite dials_data files though
    fps_copy = [
        FilePair(
            tmp_path / "tmp2.expt",
            tmp_path / "tmp2.refl",
        )
    ]
    shutil.copy(fps[0].refl, fps_copy[0].refl)
    shutil.copy(fps[0].expt, fps_copy[0].expt)
    handler.write_groupids_into_files(fps_copy)

    refls = flex.reflection_table.from_file(fps_copy[0].refl)
    assert set(refls["id"]) == {0, 1, 2, 3, 4}
    sel0 = refls["id"] == 0
    sel1 = refls["id"] == 1
    assert set(refls["group_id"].select(sel0 | sel1)) == {0}
    sel2 = refls["id"] == 2
    sel3 = refls["id"] == 3
    assert set(refls["group_id"].select(sel2 | sel3)) == {1}
    sel4 = refls["id"] == 4
    assert set(refls["group_id"].select(sel4)) == {2}

    real_example_single = f"""
---
metadata:
  timepoint:
    '{fpath}' : 1
grouping:
  group_by:
    values:
      - timepoint
    tolerances:
      - 0.1
"""
    with open(tmp_path / "real_example_single.yaml", "w") as f:
        f.write(real_example_single)

    parsed = ParsedYAML(tmp_path / "real_example_single.yaml")
    handler = get_grouping_handler(parsed, "group_by")

    fd = handler.split_files_to_groups(tmp_path, fps)
    assert list(fd.keys()) == ["group_1"]
    # with max lattices=2, 17001 has two lattices, 17002,17003,17004 have one
    filelist_1 = fd["group_1"]
    assert len(filelist_1) == 1
    expts1 = load.experiment_list(filelist_1[0].expt)
    assert len(expts1) == 5
