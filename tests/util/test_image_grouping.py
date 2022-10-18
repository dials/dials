from __future__ import annotations

import os
import subprocess
from pathlib import Path

import h5py
import numpy as np
import pytest

import dials.util.image_grouping
from dials.util.image_grouping import (  # MetadataInFile,
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

    from unittest import mock

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
    test_h5 = str(os.fspath(tmp_path)) + "/meta.h5"
    wl_example = f"""
metadata:
  timepoint:
    "/path/to/example_master.h5" : "{test_h5}:/timepoint"
  wavelength:
    "/path/to/example_master.h5" : "{test_h5}:/wavelength"
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
    tp_array = np.array([1, 0, 1, 1, 1, 1], dtype=np.int)
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

    for g in groups:
        print(g)
    ftg = GroupingImageTemplates._files_to_groups(merge_by.extract_data(), groups)
    iitgi = ftg[I1]["img_idx_to_group_id"]
    for i in range(100):
        assert iitgi[i] == i % 10
    assert iitgi.single_return_val is None
    assert iitgi.group_ids.size() == 0

    simple_block_example = """
metadata:
  crystal_id:
    "/path/to/example_#####.cbf" : "block=1:100:10"
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
    "/path/to/example_master.h5" : "/path/to/meta.h5:/timepoint"
  wavelength:
    "/path/to/example_master.h5" : "repeat=2"
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
    "/path/to/example_master.h5" : "repeat=4"
  wavelength:
    "/path/to/example_master.h5" : "repeat=2"
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


def test_real_example(tmp_path, dials_data):

    """This test tests a few use cases on real cbf data, using the template
    metadata definition.

    First, post-indexed and pre-indexed data are split (to test the different
    imageset layouts for these kinds of files), with the "repeat=" definition.
    Then, the "block=" definition is tested on the indexed data.
    Finally, the single value definition is tested.
    """

    ssx = dials_data("cunir_serial", pathlib=True)
    fpath = str(os.fspath(ssx)) + "/merlin0047_#####.cbf"
    real_example = f"""
---
metadata:
  timepoint:
    "{fpath}" : "repeat=2"
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
    groupby = parsed.groupings["group_by"]
    handler = GroupingImageTemplates(groupby)

    args = ["dials.import", f"template={fpath}"]
    result = subprocess.run(args, cwd=tmp_path, capture_output=True)
    assert not result.returncode and not result.stderr

    args = ["dials.find_spots", "imported.expt"]
    result = subprocess.run(args, cwd=tmp_path, capture_output=True)
    assert not result.returncode and not result.stderr
    args = [
        "dev.dials.ssx_index",
        "imported.expt",
        "strong.refl",
        "unit_cell=96.4,96.4,96.4,90,90,90",
        "space_group=P213",
        "nproc=1",
        "max_lattices=2",
    ]
    result = subprocess.run(args, cwd=tmp_path, capture_output=True)
    assert not result.returncode and not result.stderr
    from dxtbx.serialize import load

    # First test on indexed data - here each image has its own imageset
    # only images 17001, 17002, 17003, 17004 get indexed, so expect these to be split into groups [1,0,1,0]
    fps = [FilePair(Path(tmp_path / "indexed.expt"), Path(tmp_path / "indexed.refl"))]
    fd = handler.split_files_to_groups(tmp_path, fps)

    # with max lattices=2, 17001 has two lattices, 17002,17003,17004 have one
    assert list(fd.keys()) == ["group_1", "group_2"]
    filelist_1 = fd["group_1"]
    assert len(filelist_1) == 1
    expts1 = load.experiment_list(filelist_1[0].expt)
    for expt in expts1:
        print(expt.imageset.get_path(0).split("_")[-1])
    assert len(expts1) == 2
    assert expts1[0].imageset.get_path(0).split("_")[-1] == "17002.cbf"
    assert expts1[1].imageset.get_path(0).split("_")[-1] == "17004.cbf"
    filelist_2 = fd["group_2"]
    assert len(filelist_2) == 1
    expts2 = load.experiment_list(filelist_2[0].expt)
    assert len(expts2) == 3
    assert expts2[0].imageset.get_path(0).split("_")[-1] == "17001.cbf"
    assert expts2[1].imageset.get_path(0).split("_")[-1] == "17001.cbf"
    assert expts2[2].imageset.get_path(0).split("_")[-1] == "17003.cbf"

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
    "{fpath}" : "block=17000:17004:2"
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
    groupby = parsed.groupings["group_by"]
    handler = GroupingImageTemplates(groupby)

    fps = [FilePair(Path(tmp_path / "indexed.expt"), Path(tmp_path / "indexed.refl"))]
    fd = handler.split_files_to_groups(tmp_path, fps)
    assert list(fd.keys()) == ["group_1", "group_2", "group_3"]
    # with max lattices=2, 17001 has two lattices, 17002,17003,17004 have one
    filelist_1 = fd["group_1"]
    assert len(filelist_1) == 1
    expts1 = load.experiment_list(filelist_1[0].expt)
    assert len(expts1) == 2
    assert expts1[0].imageset.get_path(0).split("_")[-1] == "17001.cbf"
    assert expts1[1].imageset.get_path(0).split("_")[-1] == "17001.cbf"
    filelist_2 = fd["group_2"]
    assert len(filelist_2) == 1
    expts2 = load.experiment_list(filelist_2[0].expt)
    assert len(expts2) == 2
    assert expts2[0].imageset.get_path(0).split("_")[-1] == "17002.cbf"
    assert expts2[1].imageset.get_path(0).split("_")[-1] == "17003.cbf"
    filelist_3 = fd["group_3"]
    assert len(filelist_3) == 1
    expts3 = load.experiment_list(filelist_3[0].expt)
    assert len(expts3) == 1
    assert expts3[0].imageset.get_path(0).split("_")[-1] == "17004.cbf"

    handler.write_groupids_into_files(fps)
    from dials.array_family import flex

    refls = flex.reflection_table.from_file(fps[0].refl)
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
    "{fpath}" : 1
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
    groupby = parsed.groupings["group_by"]
    handler = GroupingImageTemplates(groupby)

    fd = handler.split_files_to_groups(tmp_path, fps)
    assert list(fd.keys()) == ["group_1"]
    # with max lattices=2, 17001 has two lattices, 17002,17003,17004 have one
    filelist_1 = fd["group_1"]
    assert len(filelist_1) == 1
    expts1 = load.experiment_list(filelist_1[0].expt)
    assert len(expts1) == 5
