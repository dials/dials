"""
Routines for handling arbitrary groupings of images according to metadata.

Groupings are defined in yaml files like the example yaml below, by providing images, metadata and
grouping structure.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict

import h5py
import yaml
from yaml.loader import SafeLoader

example_yaml = """
---
images:
  - "/path/to/example_master.h5"                                            # images are h5 files
  - "/path/to/example_2_master.h5"
metadata:
  timepoint:
    "/path/to/example_master.h5" : "/path/to/example_master.h5:/timepoint"  # metadata contained in image file
    "/path/to/example_2_master.h5" : "/path/to/meta.h5:/timepoint"          # metadata in separate file
  wavelength:
    /path/to/example_master.h5 : 0.4                                        # metadata is a shared value for every image
    /path/to/example_2_master.h5 : 0.6
structure:
  merge_by:                  # define a grouping for a particular process
    values:                  # the values are keys for the metadata
      - timepoint
      - wavelength
    tolerances:
      - 0.1
      - 0.01
  index_by:                  # define a grouping for a different process
    values:
      - wavelength
    tolerances:
      - 0.01
"""

simple_example = """
---
images:
  - "/path/to/example_master.h5"
metadata:
  timepoint:
    "/path/to/example_master.h5" : "/path/to/meta.h5:/timepoint"
structure:
  merge_by:
    values:
      - timepoint
    tolerances:
      - 0.1
"""

simple_template_example = """
---
templates:
  - "/path/to/example_#####.cbf"
metadata:
  dose_point:
    "/path/to/example_#####.cbf" : "repeat=10"
structure:
  merge_by:
    values:
      - dose_point
    tolerances:
      - 0.1
"""


class MetadataForFile(object):
    pass


class MetadataInFile(MetadataForFile):
    def __init__(self, file, item):
        self.file = file
        self.item = item

    def __eq__(self, other):
        return self.file == other.file and self.item == other.item

    def __str__(self):
        return f"(File={self.file}, item={self.item})"


class RepeatInImageFile(MetadataForFile):
    def __init__(self, repeat: int):
        self.repeat = repeat

    def __eq__(self, other):
        return self.repeat == other.repeat

    def __str__(self):
        return f"repeat={self.repeat}"


class ConstantMetadataForFile(MetadataForFile):
    def __init__(self, value):
        self.value = value

    def __eq__(self, other):
        return self.value == other.value

    def __str__(self):
        return str(self.value)


@dataclass(frozen=True)
class ImageFile:
    name: str
    is_h5: bool
    is_template: bool

    def __str__(self):
        if self.is_h5:
            return self.name
        else:
            return f"template={self.name}"


class ImgToMetadataDict(dict):
    # helps with type checking and printing
    def __setitem__(self, key: ImageFile, value: MetadataForFile):
        super().__setitem__(key, value)

    def __str__(self):
        return ", ".join(f"{key}: {value}" for key, value in self.items())


class NameToMetadataDict(dict):
    # helps with type checking and printing
    def __setitem__(self, key: str, value: MetadataForFile):
        super().__setitem__(key, value)

    def __str__(self):
        return ", ".join(f"{key}: {value}" for key, value in self.items())


class ParsedGrouping(object):
    def __init__(self, images: Dict[str, ImageFile], name):
        self.name = name
        self._images_to_metadata: Dict[ImageFile, NameToMetadataDict] = {
            i: NameToMetadataDict() for i in images.values()
        }
        self.tolerances: dict = {}
        self._metadata_names: set = set()

    @property
    def n_images(self) -> int:
        return len(self._images_to_metadata)

    @property
    def metadata_names(self) -> set:
        return self._metadata_names

    def add_metadata_for_image(
        self, image: ImageFile, metadata: NameToMetadataDict
    ) -> None:
        if image not in self._images_to_metadata:
            raise ValueError(f"{image} not in initialised images")
        self._images_to_metadata[image].update(metadata)
        self._metadata_names.add(list(metadata.keys())[0])

    def add_tolerances(self, tolerances: dict) -> None:
        self.tolerances = tolerances

    def check_consistent(self) -> None:
        # check we have all same keys for all images.
        if not all(self._images_to_metadata.values()):
            raise AssertionError("Metadata location only specified for some images")
        for img, v in self._images_to_metadata.items():
            if set(v.keys()) != self.metadata_names:
                raise AssertionError(
                    "Metadata names not consistent across all images:\n"
                    + f"full set of metadata names across images: {self.metadata_names}\n"
                    + f"Image {img} : metadata names: {set(v.keys())}"
                )
        if set(self.tolerances.keys()) != self.metadata_names:
            raise AssertionError(
                f"Tolerance names != metadata names: {set(self.tolerances.keys())}, {self.metadata_names}"
            )

    def __str__(self):
        tolerances = "\n".join(f"    {n} : {t}" for n, t in self.tolerances.items())
        header = f"""
Summary of data in ParsedGrouping class
  Grouping name: {self.name}
  Metadata names: {', '.join(n for n in self.metadata_names)}
  Tolerances:
{tolerances}
"""
        for i, v in self._images_to_metadata.items():
            header += f"  Image: {i}\n    metadata: {v}\n"
        return header

    def extract_data(self) -> Dict[ImageFile, Dict[str, Any]]:
        relevant_metadata: Dict[ImageFile, Dict[str, Any]] = defaultdict(dict)
        for (
            img,
            metadata_dict,
        ) in self._images_to_metadata.items():  # ImgFile, NameToMetaDict
            for k, v in metadata_dict.items():  # str, MetaForFile
                if isinstance(v, MetadataInFile):
                    # try to read the metadata from the file
                    item = v.item
                    with h5py.File(v.file, mode="r") as filedata:
                        try:
                            item = item.split("/")[1:]
                            while item:
                                next = item.pop(0)
                                filedata = filedata[next]
                            this_values = filedata[()]
                        except Exception:
                            raise ValueError(f"Unable to extract {item} from {v.file}")
                        else:
                            relevant_metadata[img][k] = this_values
                elif isinstance(v, ConstantMetadataForFile):
                    relevant_metadata[img][k] = v.value
                else:
                    raise TypeError()
        return relevant_metadata


class ParsedYAML(object):
    def __init__(self, yml_file: Path):
        with open(yml_file, "r") as f:
            data = list(yaml.load_all(f, Loader=SafeLoader))[0]

        # load the images or templates
        self._images: dict[str, ImageFile] = {}
        if ("images" not in data) and ("templates" not in data):
            raise AssertionError(
                f"No images defined in {yml_file}. Example format: {example_yaml}"
            )
        elif ("images" in data) and ("templates" in data):
            raise AssertionError(
                f"Only images or templates can be defined in {yml_file}. Example format: {example_yaml}"
            )
        elif "images" in data:
            if not isinstance(data["images"], list):
                raise AssertionError(
                    f"'images:' in {yml_file} must be defined as a list. Example format: {example_yaml}"
                )
            for name in data["images"]:
                if name.endswith(".h5") or name.endswith(".nxs"):
                    self._images[name] = ImageFile(name, True, False)
                else:
                    raise AssertionError("Image file must be .h5 or .nxs format")
        else:  # templates in data
            if not isinstance(data["templates"], list):
                raise AssertionError(
                    f"'templates:' in {yml_file} must be defined as a list. Example format: {example_yaml}"
                )
            for name in data["templates"]:
                self._images[name] = ImageFile(name, False, True)

        # Check for the metadata and structure
        if "metadata" not in data:
            raise AssertionError(
                f"No metadata defined in {yml_file}. Example format: {example_yaml}"
            )
        if not isinstance(data["metadata"], dict):
            raise AssertionError(
                f"'metadata:' in {yml_file} must be defined as a dictionary. Example format: {example_yaml}"
            )
        if "structure" not in data:
            raise AssertionError(
                f"No structure defined in {yml_file}. Example format: {example_yaml}"
            )
        if not isinstance(data["structure"], dict):
            raise AssertionError(
                f"'structure:' in {yml_file} must be defined as a dictionary. Example format: {example_yaml}"
            )
        self._yml_file = yml_file
        self.metadata_items: dict[str, ImgToMetadataDict] = {}
        # ^ e.g. timepoint to MetadataDict
        self._groupings: dict[str, ParsedGrouping] = {}
        # ^ e.g. mergeby to ParsedGrouping
        self._parse_metadata(data["metadata"])
        self._parse_structure(data["structure"])

    @property
    def groupings(self) -> Dict[str, ParsedGrouping]:
        return self._groupings

    def _parse_metadata(self, metadata: dict):
        for (
            name,
            metadict,
        ) in metadata.items():  # name is e.g. timepoint, metadict is image : file
            self.metadata_items[name] = ImgToMetadataDict()
            if not isinstance(metadict, dict):
                raise AssertionError(
                    f"Metadata items in {self._yml_file} must be defined as a dictionary of images to metadata. Example format: {example_yaml}"
                )
            for image, meta in metadict.items():
                if image not in self._images:
                    raise ValueError(
                        f"Image {image} not listed in 'images:' in {self._yml_file}"
                    )
                if type(meta) is float or type(meta) is int:
                    self.metadata_items[name][
                        self._images[image]
                    ] = ConstantMetadataForFile(meta)
                elif type(meta) is str:
                    if meta.startswith("repeat="):
                        try:
                            n = int(meta.split("=")[1])
                        except Exception as e:
                            raise ValueError(
                                f"Error interpreting {meta} as repeat=n, where n is an integer. Specific exception: {e}"
                            )
                        else:
                            self.metadata_items[name][
                                self._images[image]
                            ] = RepeatInImageFile(n)
                    else:
                        try:
                            metafile, loc = meta.split(":")
                        except Exception as e:
                            raise ValueError(
                                f"Unable to understand value: {meta}, expected format file:item e.g. /data/file.h5:/entry/data/timepoint. Specific exception: {e}"
                            )
                        else:
                            self.metadata_items[name][
                                self._images[image]
                            ] = MetadataInFile(metafile, loc)
                else:
                    raise TypeError(
                        "Only float, int and string metadata items are allowed"
                    )
        for name, items in self.metadata_items.items():
            if len(items) != len(self._images):
                raise ValueError(f"Not all images have {name} values specified")

    def _parse_structure(self, structure):
        for groupby, data in structure.items():
            self._groupings[groupby] = ParsedGrouping(self._images, groupby)
            if "values" not in data:
                raise ValueError(
                    f"Grouping {groupby} does not have 'values:' specified"
                )
            if not isinstance(data["values"], list):
                raise ValueError(
                    f"Grouping {groupby} 'values:' must be specified as a list"
                )
            if "tolerances" not in data:
                raise ValueError(
                    f"Grouping {groupby} does not have 'tolerances' specified"
                )
            if not isinstance(data["tolerances"], list):
                raise ValueError(
                    f"Grouping {groupby} 'tolerances:' must be specified as a list"
                )
            values = data["values"]
            tolerances = data["tolerances"]
            if len(tolerances) != len(values):
                raise ValueError(
                    f"The tolerances and values lists are unequal in {groupby} grouping"
                )

            for name in values:  # e.g. timepoint, wavelength
                if name not in self.metadata_items:
                    raise ValueError(
                        f"Location of {name} values not specified in 'metadata:'"
                    )
                for imagefile in self._images.values():
                    metaforname = NameToMetadataDict(
                        {name: self.metadata_items[name][imagefile]}
                    )
                    self._groupings[groupby].add_metadata_for_image(
                        imagefile, metaforname
                    )
            self._groupings[groupby].add_tolerances(
                {n: t for n, t in zip(values, tolerances)}
            )
            self._groupings[groupby].check_consistent()
