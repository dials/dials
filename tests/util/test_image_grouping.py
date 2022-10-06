from __future__ import annotations

from dials.util.image_grouping import (
    ConstantMetadataForFile,
    ImageFile,
    MetadataInFile,
    ParsedYAML,
    RepeatInImageFile,
    example_yaml,
    simple_template_example,
)


def test_yml_parsing(tmp_path):
    with open(tmp_path / "example.yaml", "w") as f:
        f.write(example_yaml)
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
        I1: MetadataInFile("/path/to/example_master.h5", "/timepoint"),
        I2: MetadataInFile("/path/to/meta.h5", "/timepoint"),
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
        "timepoint": MetadataInFile("/path/to/example_master.h5", "/timepoint"),
        "wavelength": ConstantMetadataForFile(0.4),
    }
    assert merge_by._images_to_metadata[I2] == {
        "timepoint": MetadataInFile("/path/to/meta.h5", "/timepoint"),
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


def test_yml_parsing_template(tmp_path):
    with open(tmp_path / "example.yaml", "w") as f:
        f.write(simple_template_example)
    parsed = ParsedYAML(tmp_path / "example.yaml")

    I1 = ImageFile("/path/to/example_#####.cbf", False, True)
    # check the images were found
    assert list(parsed._images.keys()) == [I1.name]

    # check the metadata was found
    assert parsed.metadata_items["dose_point"] == {I1: RepeatInImageFile(10)}

    # Now check the two groupings
    assert list(parsed.groupings.keys()) == ["merge_by"]
    merge_by = parsed.groupings["merge_by"]
    assert merge_by.metadata_names == {"dose_point"}
    assert merge_by._images_to_metadata[I1] == {"dose_point": RepeatInImageFile(10)}
    assert merge_by.tolerances == {"dose_point": 0.1}
    print(merge_by)  # test the __str__ method
