"""
Routines for handling arbitrary groupings of images according to metadata.

Groupings are defined in yaml files like the example yaml below, by providing images, metadata and
grouping structure.
"""

from __future__ import annotations

from collections import defaultdict
from copyreg import remove_extension
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, TypedDict

import h5py
import numpy as np
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
  wavelength:
    "/path/to/example_#####.cbf" : 1
structure:
  merge_by:
    values:
      - dose_point
      - wavelength
    tolerances:
      - 0.1
      - 0.01
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


@dataclass
class ExtractedValues:
    data: np.array
    is_repeat: bool


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

    def extract_data(self) -> Dict[ImageFile, Dict[str, ExtractedValues]]:
        relevant_metadata: Dict[ImageFile, Dict[str, ExtractedValues]] = defaultdict(
            dict
        )
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
                            relevant_metadata[img][k] = ExtractedValues(
                                this_values, False
                            )
                elif isinstance(v, ConstantMetadataForFile):
                    relevant_metadata[img][k] = ExtractedValues(
                        np.array([v.value]), False
                    )
                elif isinstance(v, RepeatInImageFile):
                    relevant_metadata[img][k] = ExtractedValues(
                        np.array(range(0, v.repeat)), True
                    )
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
        for (name, metadict) in metadata.items():
            # name is e.g. timepoint, metadict is image : file
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


import functools
import itertools


class MetaDataGroup(object):
    def __init__(self, data_dict):
        self._data_dict = data_dict

    def min_max_for_metadata(self, name):
        return (self._data_dict[name]["min"], self._data_dict[name]["max"])

    def __str__(self):
        return "\n".join(
            f"  {k} : {v['min']} - {v['max']}" for k, v in self._data_dict.items()
        )


class ImgIdxToGroupId(object):
    def __init__(self, single_return_val=None, repeat=None):
        self.single_return_val = single_return_val
        self.repeat = repeat
        self.group_ids = None

    def add_selection(self, int_array):
        self.group_ids = int_array

    def set_selected(self, sel, i):
        self.group_ids.set_selected(sel, i)

    def __getitem__(self, key):
        if self.single_return_val is not None:
            return self.single_return_val
        elif self.repeat is not None:
            return key % self.repeat
        return self.group_ids[key]


class GroupInfo(TypedDict):
    group_ids: List[int]
    img_idx_to_group_id: ImgIdxToGroupId


def _determine_groupings(parsed_group: ParsedGrouping):
    metadata = parsed_group.extract_data()
    tolerances = parsed_group.tolerances

    unique_values_for_metadata = {n: [] for n in parsed_group.metadata_names}
    n_images_per_file: dict[ImageFile, int] = {file: 1 for file in metadata.keys()}

    # Determine the sets of unique metadata values for each metadata name
    for name in parsed_group.metadata_names:
        values = np.array([])
        # loop over the image files to get the values of this metadata
        for file, md in metadata.items():
            extracted = md[name]
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
            extracted = md[name]
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
            sel = (full_vals >= val) & (full_vals < val + tolerances[name] + 1e-9)
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


from dxtbx import flumpy
from scitbx.array_family import flex


def _files_to_groups(
    metadata: Dict[ImageFile, Dict[str, ExtractedValues]], groups: List[MetaDataGroup]
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
            in_group = np.array([])
            # loop through metadata names, see if any data within limits
            repeat_val = None
            for n, extracted in metaforfile.items():
                data = extracted.data
                if extracted.is_repeat:
                    assert repeat_val is None
                    repeat_val = len(data)
                minv, maxv = group.min_max_for_metadata(n)
                s1 = data >= minv
                s2 = data < maxv
                if in_group.size == 0:
                    in_group = s1 & s2
                else:
                    in_group = in_group & s1 & s2
            if any(in_group):
                file_to_groups[f]["group_ids"].append(i)
                if in_group.size == 1:
                    # this means that all of the data arrays were size 1, i.e.
                    # all metadata items are simple labels, therefore all
                    # data for this image must belong to a single group.
                    file_to_groups[f]["img_idx_to_group_id"].single_return_val = i
                elif repeat_val:
                    file_to_groups[f]["img_idx_to_group_id"].repeat = repeat_val
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


class GroupsIdentifiersForExpt(object):
    def __init__(self):
        self.single_group = None
        self.groups_array = None
        self.keep_all_original = True
        self.identifiers = []
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


from dxtbx.serialize import load


def _get_expt_file_to_groupsdata(
    integrated_files: List[FilePair], img_file_to_groups: dict[ImageFile, GroupInfo]
) -> Dict[Path, GroupsIdentifiersForExpt]:

    expt_file_to_groupsdata: Dict[Path, GroupsIdentifiersForExpt] = {}

    for fp in integrated_files:
        expts = load.experiment_list(fp.expt, check_format=False)
        images = list({expt.imageset.paths()[0] for expt in expts})
        # could be all one (e.g. h5 files)
        # could be multiple within the same template - find that template perhaps?
        from dxtbx.sequence_filenames import group_files_by_imageset

        g = group_files_by_imageset(images)
        print(g)
        groupdata = GroupsIdentifiersForExpt()
        if len(g) == 1:  # the experiment list only refers to one image file/template.
            iset = list(g.keys())[0]
            image = None
            for ifile in img_file_to_groups.keys():
                if iset == ifile.name:
                    image = ifile
                    break
            if image is None:
                raise ValueError(f"Image {images[0]} not found in metadata")
            groups_for_img = img_file_to_groups[image]
            print(groups_for_img["group_ids"])
            if len(groups_for_img["group_ids"]) == 1:
                groupdata.single_group = groups_for_img["group_ids"][
                    0
                ]  # all data from this expt goes to a single group
                groupdata.unique_group_numbers = set(groups_for_img["group_ids"])
            else:
                # the image goes to several groups, we just need to know the groups
                # relevant for this subset of the image
                groups_for_this = []
                group_indices: ImgIdxToGroupId = img_file_to_groups[image][
                    "img_idx_to_group_id"
                ]
                for iset in expts.imagesets():
                    for idx in range(iset.get_array_range()):
                        groups_for_this.append(group_indices)
                    paths = iset.paths()
                    for idx in iset.indices():
                        path = paths[idx]

                for expt in expts:
                    print(dir(expt.imageset))
                    if expt.imageset.get_template():
                        print(expt.imageset.get_template())
                    print(expt.imageset.get_image_identifier())
                    img, index = expt.imageset.paths()[0], expt.imageset.indices()[0]
                    index = expt.imageset.indices()[0]
                    print(img, index)
                    idx = group_indices[index]
                    groups_for_this.append(idx)
                groupdata.groups_array = np.array(groups_for_this)
                groupdata.unique_group_numbers = set(groupdata.groups_array)
        else:
            # the expt list contains data from more than one image/template
            groups_for_this = []
            img_to_ifile = {}
            for image in images:
                for ifile in img_file_to_groups.keys():
                    if image == ifile.name:
                        img_to_ifile[image] = ifile
                        break
                if image not in img_to_ifile:
                    raise ValueError(f"Image {image} not found in metadata")
            for expt in expts:
                img, index = expt.imageset.paths()[0], expt.imageset.indices()[0]
                ifile = img_to_ifile[img]
                idx = img_file_to_groups[ifile]["img_idx_to_group_id"][index]
                groups_for_this.append(idx)
            groupdata.groups_array = np.array(groups_for_this)
            groupdata.unique_group_numbers = set(groupdata.groups_array)

        expt_file_to_groupsdata[fp.expt] = groupdata
    return expt_file_to_groupsdata
