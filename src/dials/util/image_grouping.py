"""
Routines for handling arbitrary groupings of images according to metadata.

Groupings are defined in yaml files like the example yaml below, by providing images, metadata and
grouping structure.
"""

from __future__ import annotations

import functools
import itertools
from collections import defaultdict
from dataclasses import dataclass
from multiprocessing import Pool
from pathlib import Path
from typing import Any, Dict, List, Tuple, TypedDict, Union

import h5py
import numpy as np
import yaml
from yaml.loader import SafeLoader

from dxtbx import flumpy
from dxtbx.sequence_filenames import group_files_by_imageset, template_regex
from dxtbx.serialize import load

from dials.array_family import flex

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
grouping:
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
metadata:
  timepoint:
    "/path/to/example_master.h5" : "/path/to/meta.h5:/timepoint"
grouping:
  merge_by:
    values:
      - timepoint
"""

simple_template_example = """
---
metadata:
  dose_point:
    "/path/to/example_#####.cbf" : "repeat=10"
  wavelength:
    "/path/to/example_#####.cbf" : 1
grouping:
  merge_by:
    values:
      - dose_point
      - wavelength
    tolerances:
      - 0.1
      - 0.01
"""

EPS = 1e-9

## Define classes for defining the metadata type and values/location.
# class to wrap some metadata
@dataclass
class ExtractedValues:
    data: np.array
    is_repeat: bool
    is_block: bool


class MetadataForFile(object):
    def __init__(self):
        self.extracted_data = None

    def extract(self) -> None:
        raise NotImplementedError


class MetadataInFile(MetadataForFile):
    def __init__(self, file, item):
        self.file = file
        self.item = item

    def __eq__(self, other):
        return self.file == other.file and self.item == other.item

    def __str__(self):
        return f"(File={self.file}, item={self.item})"

    def extract(self) -> None:
        item = self.item
        with h5py.File(self.file, mode="r") as filedata:
            try:
                item = item.split("/")[1:]
                while item:
                    next = item.pop(0)
                    filedata = filedata[next]
                this_values = filedata[()]
            except Exception:
                raise ValueError(f"Unable to extract {item} from {self.file}")
            else:
                self.extracted_data = ExtractedValues(this_values, False, False)


class RepeatInImageFile(MetadataForFile):
    def __init__(self, repeat: int):
        self.repeat = repeat

    def __eq__(self, other):
        return self.repeat == other.repeat

    def __str__(self):
        return f"repeat={self.repeat}"

    def extract(self) -> ExtractedValues:
        self.extracted_data = ExtractedValues(
            np.array(range(0, self.repeat)), True, False
        )


class BlockInImageFile(MetadataForFile):
    def __init__(self, first: int, last: int, block: int):
        self.first = first
        self.last = last
        self.block = block

    def __eq__(self, other: BlockInImageFile):
        return (
            (self.first == other.first)
            and (self.block == other.block)
            and (self.last == other.last)
        )

    def __str__(self):
        return f"first={self.first}, block={self.block}, last={self.last}"

    def extract(self) -> ExtractedValues:
        n_images = self.last + 1 - self.first
        n_groups = (n_images // self.block) + 1
        x = np.arange(0, stop=n_groups)
        self.extracted_data = ExtractedValues(x, False, True)
        # return ExtractedValues(np.repeat(x, self.block), False, True)


class ConstantMetadataForFile(MetadataForFile):
    def __init__(self, value):
        self.value = value

    def __eq__(self, other):
        return self.value == other.value

    def __str__(self):
        return str(self.value)

    def extract(self):
        self.extracted_data = ExtractedValues(np.array([self.value]), False, False)


@dataclass(frozen=True)  # frozen=True makes hashable to use as keys in a dict
class ImageFile:
    name: str
    is_h5: bool
    is_template: bool

    def __str__(self):
        if self.is_template:
            return f"template={self.name}"
        return self.name


class ImgToMetadataDict(dict):
    # Class to define a typed dictionary between ImageFiles and Metadata
    def __setitem__(self, key: ImageFile, value: MetadataForFile):
        super().__setitem__(key, value)

    def __str__(self):
        return ", ".join(f"{key}: {value}" for key, value in self.items())


class NameToMetadataDict(dict):
    # Class to define a typed dictionary between str (for filename) and Metadata
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
            # also don't allow mixing of range and metadata in file
            n_repeat = 0
            has_data_from_file = False
            for item in v.values():
                if isinstance(item, RepeatInImageFile):
                    n_repeat += 1
                elif isinstance(item, MetadataInFile):
                    has_data_from_file = True
            if n_repeat and has_data_from_file:
                raise AssertionError(
                    "A given image file cannot have both a repeat= definition and metadata from a file. Please add the repeat as a metadata item in the file"
                )
            elif n_repeat > 1:
                raise AssertionError(
                    "A given image file cannot have more than one repeat= definition. Please define the repeats as metadata items in a metafile"
                )
            if has_data_from_file and img.is_template:
                raise AssertionError(
                    "Template definitions can only be paired with metadata of the form 'repeat=' or single values per template"
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

    def extract_data(self) -> Dict[ImageFile, Dict[str, MetadataForFile]]:
        relevant_metadata: Dict[ImageFile, Dict[str, MetadataForFile]] = defaultdict(
            dict
        )
        for img, meta in self._images_to_metadata.items():  # ImgFile, NameToMetaDict
            for name, metaforfile in meta.items():  # str, MetaForFile
                metaforfile.extract()
                relevant_metadata[img][name] = metaforfile
        return relevant_metadata


class ParsedYAML(object):
    def __init__(self, yml_file: Path):
        with open(yml_file, "r") as f:
            data = list(yaml.load_all(f, Loader=SafeLoader))[0]
        # load the images or templates
        self._images: dict[str, ImageFile] = {}

        # Check for the metadata first
        if "metadata" not in data:
            raise AssertionError(
                f"No metadata defined in {yml_file}. Example format: {example_yaml}"
            )
        if not isinstance(data["metadata"], dict):
            raise AssertionError(
                f"'metadata:' in {yml_file} must be defined as a dictionary. Example format: {example_yaml}"
            )
        if "images" in data:
            if not isinstance(data["images"], list):
                raise AssertionError(
                    f"'images:' in {yml_file} must be defined as a list. Example format: {example_yaml}"
                )
            for name in data["images"]:
                if name.endswith(".h5") or name.endswith(".nxs"):
                    self._images[name] = ImageFile(name, True, False)
                else:
                    raise AssertionError("Image file must be .h5 or .nxs format")
        elif "templates" in data:
            if not isinstance(data["templates"], list):
                raise AssertionError(
                    f"'templates:' in {yml_file} must be defined as a list. Example format: {example_yaml}"
                )
            for name in data["templates"]:
                self._images[name] = ImageFile(name, False, True)
        else:
            self._images = self._extract_images_from_metadata(data["metadata"])

        if "grouping" not in data:
            raise AssertionError(
                f"No grouping defined in {yml_file}. Example format: {example_yaml}"
            )
        if not isinstance(data["grouping"], dict):
            raise AssertionError(
                f"'grouping:' in {yml_file} must be defined as a dictionary. Example format: {example_yaml}"
            )
        self._yml_file = yml_file
        self.metadata_items: dict[str, ImgToMetadataDict] = {}
        # ^ e.g. timepoint to MetadataDict
        self._groupings: dict[str, ParsedGrouping] = {}
        # ^ e.g. mergeby to ParsedGrouping
        self._parse_metadata(data["metadata"])
        self._parse_grouping_structure(data["grouping"])

    def _extract_images_from_metadata(self, metadata: dict) -> dict[str, ImageFile]:
        images: dict[str, ImageFile] = {}
        for i, metadict in enumerate(metadata.values()):
            # name is e.g. timepoint, metadict is image : file
            if not isinstance(metadict, dict):
                raise ValueError(
                    f"Metadata items in {self._yml_file} must be defined as a dictionary of images to metadata. Example format: {example_yaml}"
                )
            if i == 0:
                for image in metadict.keys():
                    if "#" in image:
                        images[image] = ImageFile(image, False, True)
                    elif image.endswith(".h5") or image.endswith(".nxs"):
                        images[image] = ImageFile(image, True, False)
                    else:
                        raise ValueError(
                            "Image file must be .h5 or .nxs format, or be an image template (containing #)"
                        )
                if any(i.is_h5 for i in images.values()) and any(
                    i.is_template for i in images.values()
                ):
                    raise ValueError(
                        "Cannot mix image file and templates as definitions of images"
                    )
            else:
                new_images: dict[str, ImageFile] = {}
                for image in metadict.keys():
                    if "#" in image:
                        new_images[image] = ImageFile(image, False, True)
                    elif image.endswith(".h5") or image.endswith(".nxs"):
                        new_images[image] = ImageFile(image, True, False)
                    else:
                        raise ValueError(
                            "Image file must be .h5 or .nxs format, or be an image template (containing #)"
                        )
                    if new_images != images:
                        raise ValueError(
                            "Inconsistent images for different metadata items"
                        )
        return images

    @property
    def groupings(self) -> Dict[str, ParsedGrouping]:
        return self._groupings

    def _parse_metadata(self, metadata: dict):
        for (name, metadict) in metadata.items():
            # name is e.g. timepoint, metadict is image : file
            self.metadata_items[name] = ImgToMetadataDict()
            if not isinstance(metadict, dict):
                raise AssertionError(
                    f"Metadata items in {self._yml_file} must be defined as a dictionary of images to metadata. Example format: {example_yaml}"
                )
            for image, meta in metadict.items():
                try:
                    imgfile = self._images[image]
                except KeyError:
                    raise ValueError(
                        f"Image {image} not listed in 'images:' in {self._yml_file}"
                    )
                if type(meta) is float or type(meta) is int:
                    self.metadata_items[name][imgfile] = ConstantMetadataForFile(meta)
                elif type(meta) is str:
                    if meta.startswith("repeat="):
                        try:
                            n = int(meta.split("=")[1])
                        except Exception as e:
                            raise ValueError(
                                f"Error interpreting {meta} as repeat=n, where n is an integer. Specific exception: {e}"
                            )
                        else:
                            self.metadata_items[name][imgfile] = RepeatInImageFile(n)
                    elif meta.startswith("block="):
                        try:
                            first, last, block = meta.split("=")[1].split(":")
                            first = int(first)
                            last = int(last)
                            block = int(block)
                        except Exception as e:
                            raise ValueError(
                                f"Error interpreting {meta} as block=first:last:n, where first,last are the first and last image numbers and n is an"
                                + f"\ninteger indicating the number of images in a block. Specific exception: {e}"
                            )
                        else:
                            self.metadata_items[name][imgfile] = BlockInImageFile(
                                first, last, block
                            )
                    else:
                        try:
                            metafile, loc = meta.split(":")
                        except Exception as e:
                            raise ValueError(
                                f"Unable to understand value: {meta}, expected format file:item e.g. /data/file.h5:/entry/data/timepoint. Specific exception: {e}"
                            )
                        else:
                            self.metadata_items[name][imgfile] = MetadataInFile(
                                metafile, loc
                            )
                else:
                    raise TypeError(
                        "Only float, int and string metadata items are allowed"
                    )

        for name, items in self.metadata_items.items():
            if len(items) != len(self._images):
                raise ValueError(f"Not all images have {name} values specified")

    def _parse_grouping_structure(self, structure):
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
            values = data["values"]
            if "tolerances" in data:
                if not isinstance(data["tolerances"], list):
                    raise ValueError(
                        f"Grouping {groupby} 'tolerances:' must be specified as a list"
                    )
                tolerances = data["tolerances"]
                if len(tolerances) != len(values):
                    raise ValueError(
                        f"The tolerances and values lists are unequal in {groupby} grouping: {tolerances}, {values}"
                    )
            else:
                tolerances = [0 for _ in values]

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


# Class to store metadata group parameters
class MetaDataGroup(object):
    def __init__(self, data_dict):
        self._data_dict: dict[str, dict[str, float]] = data_dict
        self._default_all = False

    def min_max_for_metadata(self, name):
        return (self._data_dict[name]["min"], self._data_dict[name]["max"])

    def __str__(self):
        if self._default_all:
            return "all data"
        outlines = []
        for k, v in self._data_dict.items():
            if abs(v["min"] - v["max"]) < EPS:
                outlines.append(f"{k}={v['min']}")
            else:
                outlines.append(f"{k}={v['min']}-{v['max']}")
        return ", ".join(outlines)


# Define mapping from image index to group id.
class ImgIdxToGroupId(object):
    def __init__(self, single_return_val=None, repeat=None, block=None):
        self.single_return_val = single_return_val
        self.repeat = repeat
        self.group_ids = None  # only relevant for h5 files - split up?
        self.block = None

    def add_selection(self, int_array):
        self.group_ids = int_array

    def set_selected(self, sel, i):
        self.group_ids.set_selected(sel, i)

    def __getitem__(self, key):
        if self.single_return_val is not None:
            return self.single_return_val
        elif self.repeat is not None:
            return key % self.repeat
        elif self.block is not None:
            return (key - self.block.first) // self.block.block
        return self.group_ids[key]


class GroupInfo(TypedDict):
    group_ids: List[int]
    img_idx_to_group_id: ImgIdxToGroupId


def _determine_groupings(parsed_group: ParsedGrouping):
    """Determine the unique groups into which the data should be split."""
    metadata = parsed_group.extract_data()
    tolerances = parsed_group.tolerances

    unique_values_for_metadata = {n: [] for n in parsed_group.metadata_names}
    n_images_per_file: dict[ImageFile, int] = {file: 1 for file in metadata.keys()}

    # Determine the sets of unique metadata values for each metadata name
    for name in parsed_group.metadata_names:
        values = np.array([])
        # loop over the image files to get the values of this metadata
        for file, md in metadata.items():
            extracted = md[name].extracted_data
            values = np.concatenate([values, extracted.data])
            n_images_per_file[file] = extracted.data.size
        set_of_values = np.array(sorted(set(np.around(values, decimals=9))))
        if len(set_of_values) == 1:
            unique_values_for_metadata[name] = set_of_values
        else:
            diffs = np.abs(set_of_values[1:] - set_of_values[:-1])
            if np.min(diffs) > tolerances[name]:
                unique_values_for_metadata[name] = set_of_values
            else:
                unique_vals = [set_of_values[0]]
                for tp in set_of_values[1:]:
                    if tp - unique_vals[-1] > tolerances[name]:
                        unique_vals.append(tp)
                unique_values_for_metadata[name] = np.array(unique_vals)

    groups: List[MetaDataGroup] = []  # list of dicts
    # NB this is not necessarily the number of actual images, in the case
    # that the metadata for a file is only single values rather than arrays.
    n_images = sum(n_images_per_file.values())
    # Create a single array of values for the effective number of images, needed
    # for doing repeated selections below when determining group membership.
    full_values_per_metadata: dict[str, np.array] = {}
    for name in parsed_group.metadata_names:
        values = np.array([])
        for file, md in metadata.items():
            extracted = md[name].extracted_data
            values = np.concatenate([values, extracted.data])
        full_values_per_metadata[name] = values

    # Now work out the all combinations of the metadata groups, and check which
    # of these are actually populated. Record the valid ranges for each group.
    combs = list(itertools.product(*unique_values_for_metadata.values()))
    for vals in combs:
        sel1 = np.full((n_images,), True)
        for name, val in zip(
            parsed_group.metadata_names, vals
        ):  # val is the lower bound for that group
            full_vals = full_values_per_metadata[name]
            sel = (full_vals >= val) & (full_vals < val + tolerances[name] + EPS)
            sel1 = sel1 & sel
        if np.any(sel1):
            groups.append(
                MetaDataGroup(
                    {
                        n: {"min": v, "max": v + tolerances[n]}
                        for n, v in zip(parsed_group.metadata_names, vals)
                    }
                )
            )
    return groups


class GroupsIdentifiersForExpt(object):
    def __init__(self):
        self.single_group = None
        self.groups_array = None
        self.unique_group_numbers = None


@dataclass(eq=False)
class FilePair:
    expt: Path
    refl: Path

    def check(self):
        if not self.expt.is_file():
            raise FileNotFoundError(f"File {self.expt} does not exist")
        if not self.refl.is_file():
            raise FileNotFoundError(f"File {self.refl} does not exist")

    def validate(self):
        expt = load.experiment_list(self.expt, check_format=False)
        refls = flex.reflection_table.from_file(self.refl)
        refls.assert_experiment_identifiers_are_consistent(expt)

    def __eq__(self, other):
        if self.expt == other.expt and self.refl == other.refl:
            return True
        return False


@dataclass
class InputIterable(object):
    working_directory: Path
    fp: FilePair
    fileindex: int
    groupindex: int
    groupdata: GroupsIdentifiersForExpt
    name: str
    reduction_params: Any


def save_subset(input_: InputIterable) -> Union[Tuple[str, FilePair], None]:
    expts = load.experiment_list(input_.fp.expt, check_format=False)
    refls = flex.reflection_table.from_file(input_.fp.refl)
    groupdata = input_.groupdata
    if (groupdata.single_group is not None) and (
        groupdata.single_group == input_.groupindex
    ):
        pass
    else:
        # need to select
        identifiers = expts.identifiers()
        sel = input_.groupdata.groups_array == input_.groupindex
        sel_identifiers = list(identifiers.select(flumpy.from_numpy(sel)))
        expts.select_on_experiment_identifiers(sel_identifiers)
        refls = refls.select_on_experiment_identifiers(sel_identifiers)
    if expts:
        exptout = (
            input_.working_directory
            / f"group_{input_.groupindex}_{input_.fileindex}.expt"
        )
        reflout = (
            input_.working_directory
            / f"group_{input_.groupindex}_{input_.fileindex}.refl"
        )
        expts.as_file(exptout)
        refls.as_file(reflout)
        return (input_.name, FilePair(exptout, reflout))
    return None


class GroupingImageTemplates(object):

    """Class that takes a parsed group and determines the groupings and mappings
    required to split input data into groups.

    This class provides specific implementations for when the images are provided
    as a template. The main difference from h5 images is getting the image index.
    """

    def __init__(self, parsed_group: ParsedGrouping, nproc=1):
        self._parsed_group = parsed_group
        self.nproc = nproc
        self._grouping_metadata = parsed_group.extract_data()
        self._groups = _determine_groupings(parsed_group)
        self._files_to_groups_dict = self._files_to_groups(
            self._grouping_metadata, self._groups
        )

    @staticmethod
    def _files_to_groups(
        metadata: Dict[ImageFile, Dict[str, MetadataForFile]],
        groups: List[MetaDataGroup],
    ) -> dict[ImageFile, GroupInfo]:

        # Ok now we have the groupings of the metadata. Now find which groups each
        # file contains.
        # Purpose here is to create an object that will allow easy allocation from
        # image to group
        file_to_groups: dict[ImageFile, GroupInfo] = {
            n: {"group_ids": [], "img_idx_to_group_id": ImgIdxToGroupId()}
            for n in metadata.keys()
        }
        for f in file_to_groups:
            metaforfile = metadata[f]
            for i, group in enumerate(groups):
                in_group = np.array([True], dtype=bool)
                # loop through metadata names, see if any data within limits
                block = None
                repeat = None
                special_types = None
                for n, extracted in metaforfile.items():  # e.g. timepoint, wavelength
                    extractedvalues = extracted.extracted_data
                    if extractedvalues.is_repeat or extractedvalues.is_block:
                        if special_types:  # i.e. already have one of these cases
                            raise ValueError(
                                "Templates only support one instance of repeat metadata or block per grouping"
                            )
                        if extractedvalues.is_repeat:
                            special_types = "repeat"
                            repeat = len(extractedvalues.data)
                        else:
                            special_types = "block"
                            block = extracted
                    # allowable data types are block, repeat or single valiue
                    minv, maxv = group.min_max_for_metadata(n)
                    s1 = extractedvalues.data >= minv
                    s2 = extractedvalues.data < maxv + EPS
                    in_group = (s1 & s2) & in_group
                    if not any(in_group):
                        break
                if any(in_group):
                    file_to_groups[f]["group_ids"].append(i)
                    if not special_types:
                        # this means that all of the data arrays were size 1, i.e.
                        # all metadata items are simple labels, therefore all
                        # data for this image must belong to a single group.
                        file_to_groups[f]["img_idx_to_group_id"].single_return_val = i
                    elif special_types == "repeat":
                        file_to_groups[f]["img_idx_to_group_id"].repeat = repeat
                    else:
                        file_to_groups[f]["img_idx_to_group_id"].block = block

        return file_to_groups

    def get_expt_file_to_groupsdata(self, data_file_pairs: List[FilePair]):
        expt_file_to_groupsdata: Dict[Path, GroupsIdentifiersForExpt] = {}

        for fp in data_file_pairs:
            expts = load.experiment_list(fp.expt, check_format=False)
            # need to match the images to the imagesets.
            images = set()
            for iset in expts.imagesets():
                images.update(iset.paths())

            template_to_group_indices = {}
            for iset in group_files_by_imageset(images):
                image = None
                for ifile in self._files_to_groups_dict.keys():
                    if iset == ifile.name:
                        image = ifile
                        break
                if image is None:
                    raise ValueError(f"Imageset {iset} not found in metadata")
                template_to_group_indices[iset] = self._files_to_groups_dict[image][
                    "img_idx_to_group_id"
                ]

            groupdata = GroupsIdentifiersForExpt()

            if (
                len(template_to_group_indices) == 1
            ):  # the experiment list only refers to one template.
                group_indices: ImgIdxToGroupId = list(
                    template_to_group_indices.values()
                )[0]
                if group_indices.single_return_val is not None:
                    groupdata.single_group = group_indices.single_return_val
                    groupdata.unique_group_numbers = {group_indices.single_return_val}
                else:
                    # the image goes to several groups, we just need to know the groups
                    # relevant for these images
                    groups_for_this = []
                    for expt in expts:
                        if expt.scan:
                            start = expt.scan.get_image_range()[0]
                            groups_for_this.append(group_indices[start])
                        else:
                            p = expt.imageset.paths()[0]
                            t = template_regex(p)
                            groups_for_this.append(group_indices[t[1]])
                    groupdata.groups_array = np.array(groups_for_this, dtype=np.uint64)
                    groupdata.unique_group_numbers = set(groupdata.groups_array)
            else:
                # the expt list contains data from more than one image/template
                groups_for_this = []
                for expt in expts:
                    if expt.scan:
                        p = expt.imageset.paths()[0]
                        t = template_regex(p)
                        group_indices = template_to_group_indices[t[0]]
                        start = expt.scan.get_image_range()[0]
                        groups_for_this.append(group_indices[start])
                    else:
                        p = expt.imageset.paths()[0]
                        t = template_regex(p)
                        group_indices = template_to_group_indices[t[0]]
                        groups_for_this.append(group_indices[t[1]])
                groupdata.groups_array = np.array(groups_for_this, dtype=np.uint64)
                groupdata.unique_group_numbers = set(groupdata.groups_array)
            expt_file_to_groupsdata[fp.expt] = groupdata
        return expt_file_to_groupsdata

    def write_groupids_into_files(self, data_file_pairs: List[FilePair]) -> None:
        "Write a group_id column into the reflection table"
        expt_file_to_groupsdata: Dict[
            Path, GroupsIdentifiersForExpt
        ] = self.get_expt_file_to_groupsdata(data_file_pairs)

        def set_group_id_column(
            filepair: FilePair,
            groupdata: GroupsIdentifiersForExpt,
        ):
            expts = load.experiment_list(filepair.expt, check_format=False)
            refls = flex.reflection_table.from_file(filepair.refl)
            groupdata = groupdata
            if groupdata.single_group is not None:
                group_id = groupdata.single_group
                refls["group_id"] = flex.int(refls.size(), group_id)
            else:
                identifiers_map = refls.experiment_identifiers()
                list_of_identifiers: List = list(expts.identifiers())
                groups_array = flumpy.from_numpy(groupdata.groups_array)

                def table_id_to_group_id(table_id):
                    identifier = identifiers_map[table_id]
                    idx = list_of_identifiers.index(identifier)
                    return groups_array[idx]

                refls["group_id"] = flex.int(map(table_id_to_group_id, refls["id"]))

            refls.as_file(filepair.refl)
            return filepair

        for fp in data_file_pairs:
            groupdata = expt_file_to_groupsdata[fp.expt]
            set_group_id_column(fp, groupdata)

    def split_files_to_groups(
        self,
        working_directory: Path,
        data_file_pairs: List[FilePair],
        reduction_params,
        function_to_apply=save_subset,
        prefix="",
    ):

        expt_file_to_groupsdata: Dict[
            Path, GroupsIdentifiersForExpt
        ] = self.get_expt_file_to_groupsdata(data_file_pairs)
        template = "{name}group_{index:0{maxindexlength:d}d}"
        name_template = functools.partial(
            template.format,
            name=f"{prefix}",
            maxindexlength=len(str(len(self._groups))),
        )

        names: List[str] = [
            name_template(index=i + 1) for i, _ in enumerate(self._groups)
        ]
        filesdict: dict[str, List[FilePair]] = {name: [] for name in names}

        input_iterable = []
        for groupindex, name in enumerate(names):
            for fileindex, fp in enumerate(data_file_pairs):
                groupdata = expt_file_to_groupsdata[fp.expt]
                if groupindex in groupdata.unique_group_numbers:
                    input_iterable.append(
                        InputIterable(
                            working_directory,
                            fp,
                            fileindex,
                            groupindex,
                            groupdata,
                            name,
                            reduction_params,
                        )
                    )
        if input_iterable:
            with Pool(min(self.nproc, len(input_iterable))) as pool:
                results = pool.map(function_to_apply, input_iterable)
            for result in results:
                if result:
                    name = result[0]
                    fp = result[1]
                    filesdict[name].append(fp)
        return filesdict


class GroupingImageFiles(GroupingImageTemplates):

    """This class provides specific implementations for when the images are h5 files.
    The main difference from tempaltes is getting the image index.
    """

    @staticmethod
    def _files_to_groups(
        metadata: Dict[ImageFile, Dict[str, MetadataForFile]],
        groups: List[MetaDataGroup],
    ) -> dict[ImageFile, GroupInfo]:

        # Ok now we have the groupings of the metadata. Now find which groups each
        # file contains.
        # Purpose here is to create an object that will allow easy allocation from
        # image to group
        file_to_groups: dict[ImageFile, GroupInfo] = {
            n: {"group_ids": [], "img_idx_to_group_id": ImgIdxToGroupId()}
            for n in metadata.keys()
        }
        for f in file_to_groups:
            metaforfile = metadata[f]
            for i, group in enumerate(groups):
                # Do we know how many images are in this h5 file?
                # N.B. assumes that data array was written with correct size.
                known_size = None
                for v in metaforfile.values():
                    if isinstance(v, MetadataInFile):
                        if v.extracted_data is not None:
                            known_size = v.extracted_data.data.size

                in_group = np.array([True], dtype=bool)  # important that this is size 1
                # loop through metadata names, see if any data within limits
                block = None
                repeat = None
                special_types = None
                for n, extracted in metaforfile.items():  # e.g. timepoint, wavelength
                    extractedvalues = extracted.extracted_data
                    if extractedvalues.is_repeat or extractedvalues.is_block:
                        if (
                            special_types or known_size
                        ):  # i.e. already have one of these cases
                            raise ValueError(
                                "H5 images only support one instance of repeat metadata or block per grouping"
                            )
                        if extractedvalues.is_repeat:
                            special_types = "repeat"
                            repeat = len(extractedvalues.data)
                        else:
                            special_types = "block"
                            block = extracted
                    minv, maxv = group.min_max_for_metadata(n)
                    s1 = extractedvalues.data >= minv
                    s2 = extractedvalues.data < maxv + EPS
                    # N.B. The checks above ensure that if s1 & s2 are greater
                    # than size 1 (i.e. size n), then in future iterations of the
                    # loop they will also be either size 1 or n, thus allowing
                    # multiplication as numpy arrays.
                    in_group = (s1 & s2) & in_group
                    if not any(in_group):
                        break
                if any(in_group):
                    file_to_groups[f]["group_ids"].append(i)
                    if not special_types and not known_size:
                        # this means that all of the data arrays were size 1, i.e.
                        # all metadata items are simple labels, therefore all
                        # data for this image must belong to a single group.
                        file_to_groups[f]["img_idx_to_group_id"].single_return_val = i
                    elif special_types == "repeat":
                        file_to_groups[f]["img_idx_to_group_id"].repeat = repeat
                    elif special_types == "block":
                        file_to_groups[f]["img_idx_to_group_id"].block = block
                    else:
                        if file_to_groups[f]["img_idx_to_group_id"].group_ids:
                            file_to_groups[f]["img_idx_to_group_id"].set_selected(
                                flumpy.from_numpy(in_group), i
                            )
                        else:
                            file_to_groups[f]["img_idx_to_group_id"].add_selection(
                                flex.int(in_group.size, 0)
                            )
                            file_to_groups[f]["img_idx_to_group_id"].set_selected(
                                flumpy.from_numpy(in_group), i
                            )

        return file_to_groups

    def get_expt_file_to_groupsdata(self, data_file_pairs: List[FilePair]):
        expt_file_to_groupsdata: Dict[Path, GroupsIdentifiersForExpt] = {}

        for fp in data_file_pairs:
            expts = load.experiment_list(fp.expt, check_format=False)
            # need to match the images to the imagesets.
            images = set()
            image_to_group_info = {}
            for iset in expts.imagesets():
                images.update(iset.paths())
            for iset in images:
                image = None
                for ifile in self._files_to_groups_dict.keys():
                    if iset == ifile.name:
                        image = ifile
                        break
                if image is None:
                    raise ValueError(f"Imageset {iset} not found in metadata")
                image_to_group_info[iset] = self._files_to_groups_dict[image]

            groupdata = GroupsIdentifiersForExpt()

            if len(images) == 1:
                group_info: GroupInfo = image_to_group_info[list(images)[0]]
                if len(group_info["group_ids"]) == 1:
                    groupdata.single_group = group_info["group_ids"][
                        0
                    ]  # all data from this expt goes to a single group
                    groupdata.unique_group_numbers = set(group_info["group_ids"])
                else:  # one h5 image, but more than one group
                    groups_for_this = []
                    group_indices = group_info["img_idx_to_group_id"]
                    for expt in expts:
                        index = expt.imageset.indices()[0]
                        idx = group_indices[index]
                        groups_for_this.append(idx)
                    groupdata.groups_array = np.array(groups_for_this, dtype=np.uint64)
                    groupdata.unique_group_numbers = set(groupdata.groups_array)
            else:  # multiple h5 images
                groups_for_this = []
                for expt in expts:
                    img, index = expt.imageset.paths()[0], expt.imageset.indices()[0]
                    idx = image_to_group_info[img]["img_idx_to_group_id"][index]
                    groups_for_this.append(idx)
                groupdata.groups_array = np.array(groups_for_this, dtype=np.uint64)
                groupdata.unique_group_numbers = set(groupdata.groups_array)

            expt_file_to_groupsdata[fp.expt] = groupdata
        return expt_file_to_groupsdata
