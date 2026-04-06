from __future__ import annotations

import os
import re

import libtbx.phil
from dxtbx.model.experiment_list import ExperimentListFactory
from libtbx.utils import Sorry

from dials.array_family import flex

_UNLOADED = object()  # sentinel for "data not yet loaded"


class FilenameDataWrapper:
    """
    Wraps a filename and its lazily-loaded data object.

    Data is loaded on first access of the `.data` property so that
    constructing a wrapper (e.g. during diff_phil string comparison) does
    not trigger a disk read.  Bad paths are detected eagerly at construction
    time via os.stat so user-facing errors still surface at argument-parse
    time.
    """

    __slots__ = ("filename", "_data", "_loader")

    def __init__(self, filename, data=_UNLOADED, _loader=None):
        self.filename = filename
        # Allow callers to pass pre-loaded data (e.g. from image-file import).
        self._data = data
        # _loader is a zero-argument callable that returns the data object.
        self._loader = _loader

    @property
    def data(self):
        if self._data is _UNLOADED:
            if self._loader is None:
                raise RuntimeError(f"No loader available for {self.filename!r}")
            self._data = self._loader()
        return self._data

    @data.setter
    def data(self, value):
        self._data = value

    def __iter__(self):
        """Support two-tuple unpacking: filename, data = wrapper."""
        yield self.filename
        yield self.data

    def __repr__(self):
        return f"FilenameDataWrapper(filename={self.filename!r})"


class ExperimentListConverters:
    """A phil converter for the experiment list class."""

    phil_type = "experiment_list"

    def __init__(self, check_format=True):
        self._check_format = check_format

    def __str__(self):
        return self.phil_type

    def from_words(self, words, master):
        s = libtbx.phil.str_from_words(words=words)
        if s is None:
            return None
        if s == "<image files>":
            return FilenameDataWrapper(filename=s, data=None)
        if not os.path.exists(s):
            raise Sorry(f"File {s} does not exist")
        # Capture check_format now; load lazily on first .data access.
        check_format = self._check_format
        abspath = os.path.abspath(s)
        return FilenameDataWrapper(
            filename=s,
            _loader=lambda: ExperimentListFactory.from_json_file(
                abspath, check_format=check_format
            ),
        )

    def as_words(self, python_object, master):
        if python_object is None:
            value = "None"
        else:
            value = python_object.filename
        return [libtbx.phil.tokenizer.word(value=value)]


class ReflectionTableConverters:
    """A phil converter for the reflection table class."""

    phil_type = "reflection_table"

    def __str__(self):
        return self.phil_type

    def from_words(self, words, master):
        s = libtbx.phil.str_from_words(words=words)
        if s is None:
            return None
        if not os.path.exists(s):
            raise Sorry(f"File {s} does not exist")
        # Load lazily on first .data access.
        abspath = os.path.abspath(s)
        return FilenameDataWrapper(
            filename=s,
            _loader=lambda: flex.reflection_table.from_file(abspath),
        )

    def as_words(self, python_object, master):
        if python_object is None:
            value = "None"
        else:
            value = python_object.filename
        return [libtbx.phil.tokenizer.word(value=value)]


class ReflectionTableSelectorConverters:
    """A phil converter for the reflection table selector class."""

    phil_type = "reflection_table_selector"

    def __str__(self):
        return self.phil_type

    def from_words(self, words, master):
        s = libtbx.phil.str_from_words(words=words)
        if s is None:
            return None
        regex = r"^\s*([\w\.]+)\s*(<=|!=|==|>=|<|>|&)\s*(.+)\s*$"
        matches = re.findall(regex, s)
        if len(matches) == 0:
            raise RuntimeError(f"{s} is not a valid selection")
        elif len(matches) > 1:
            raise RuntimeError(f"{s} is not a single selection")
        else:
            match = matches[0]
        assert len(match) == 3
        col, op, value = match
        return flex.reflection_table_selector(col, op, value)

    def as_words(self, python_object, master):
        if python_object is None:
            value = "None"
        else:
            value = (
                f"{python_object.column}{python_object.op_string}{python_object.value}"
            )
        return [libtbx.phil.tokenizer.word(value=value)]


# Create the default converter registry with the extract converters
default_converter_registry = libtbx.phil.extended_converter_registry(
    additional_converters=[
        ExperimentListConverters,
        ReflectionTableConverters,
        ReflectionTableSelectorConverters,
    ]
)


def parse(
    input_string=None,
    source_info=None,
    file_name=None,
    converter_registry=None,
    process_includes=False,
):
    """Redefinition of the parse function."""
    if converter_registry is None:
        converter_registry = default_converter_registry
    return libtbx.phil.parse(
        input_string=input_string,
        source_info=source_info,
        file_name=file_name,
        converter_registry=converter_registry,
        process_includes=process_includes,
    )
