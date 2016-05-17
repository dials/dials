#!/usr/bin/env python
#
# phil.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division
import libtbx.phil


class FilenameDataWrapper(object):
  ''' A wrapper class to store data with a filename. '''

  def __init__(self, filename, data):
    self.filename = filename
    self.data = data


class DataBlockConverters(object):
  ''' A phil converter for datablocks. '''

  phil_type = "datablock"

  cache = {}

  def __init__(self, check_format=True):
    self._check_format = check_format

  def __str__(self):
    return self.phil_type

  def from_string(self, s):
    from dxtbx.datablock import DataBlockFactory
    from os.path import exists
    from libtbx.utils import Sorry
    if s is None:
      return None
    if s not in self.cache:
      if not exists(s):
        raise Sorry('File %s does not exist' % s)
      self.cache[s] = FilenameDataWrapper(s,
        DataBlockFactory.from_json_file(s,
          check_format=self._check_format))
    return self.cache[s]

  def from_words(self, words, master):
    return self.from_string(libtbx.phil.str_from_words(words=words))

  def as_words(self, python_object, master):
    if python_object is None:
      value = "None"
    else:
      value = python_object.filename
    return [libtbx.phil.tokenizer.word(value=value)]


class ExperimentListConverters(object):
  ''' A phil converter for the experiment list class. '''

  phil_type = "experiment_list"

  cache = {}

  def __init__(self, check_format=True):
    self._check_format = check_format

  def __str__(self):
    return self.phil_type

  def from_string(self, s):
    from dxtbx.model.experiment.experiment_list import ExperimentListFactory
    from os.path import exists
    from libtbx.utils import Sorry
    if s is None:
      return None
    if s not in self.cache:
      if not exists(s):
        raise Sorry('File %s does not exist' % s)
      self.cache[s] = FilenameDataWrapper(s,
        ExperimentListFactory.from_json_file(s,
          check_format=self._check_format))
    return self.cache[s]

  def from_words(self, words, master):
    return self.from_string(libtbx.phil.str_from_words(words=words))

  def as_words(self, python_object, master):
    if python_object is None:
      value = "None"
    else:
      value = python_object.filename
    return [libtbx.phil.tokenizer.word(value=value)]


class ReflectionTableConverters(object):
  ''' A phil converter for the reflection table class. '''

  phil_type = "reflection_table"

  cache = {}

  def __str__(self):
    return self.phil_type

  def from_string(self, s):
    from dials.array_family import flex
    from os.path import exists
    from libtbx.utils import Sorry
    if s is None:
      return None
    if s not in self.cache:
      if not exists(s):
        raise Sorry('File %s does not exist' % s)
      self.cache[s] = FilenameDataWrapper(s, flex.reflection_table.from_pickle(s))
    return self.cache[s]

  def from_words(self, words, master):
    return self.from_string(libtbx.phil.str_from_words(words=words))

  def as_words(self, python_object, master):
    if python_object is None:
      value = "None"
    else:
      value = python_object.filename
    return [libtbx.phil.tokenizer.word(value=value)]


class ReflectionTableSelectorConverters(object):
  ''' A phil converter for the reflection table selector class. '''

  phil_type = "reflection_table_selector"

  cache = {}

  def __str__(self):
    return self.phil_type

  def from_string(self, s):
    from dials.array_family import flex
    if s is None:
      return None
    import re
    regex = r"^\s*([\w\.]+)\s*(<=|!=|==|>=|<|>|&)\s*(.+)\s*$"
    matches = re.findall(regex, s)
    if len(matches) == 0:
      raise RuntimeError('%s is not a valid selection' % s)
    elif len(matches) > 1:
      raise RuntimeError('%s is not a single selection' % s)
    else:
      match = matches[0]
    assert(len(match) == 3)
    col, op, value = match
    return flex.reflection_table_selector(col, op, value)

  def from_words(self, words, master):
    return self.from_string(libtbx.phil.str_from_words(words=words))

  def as_words(self, python_object, master):
    if python_object is None:
      value = "None"
    else:
      value = "%s%s%s" % (python_object.column,
                          python_object.op_string,
                          python_object.value)
    return [libtbx.phil.tokenizer.word(value=value)]


# Create the default converter registry with the extract converters
default_converter_registry = libtbx.phil.extended_converter_registry(
  additional_converters=[
    DataBlockConverters,
    ExperimentListConverters,
    ReflectionTableConverters,
    ReflectionTableSelectorConverters
  ])


def parse(
      input_string=None,
      source_info=None,
      file_name=None,
      converter_registry=None,
      process_includes=False):
  ''' Redefinition of the parse function. '''
  if converter_registry is None:
    converter_registry = default_converter_registry
  return libtbx.phil.parse(
    input_string=input_string,
    source_info=source_info,
    file_name=file_name,
    converter_registry=converter_registry,
    process_includes=process_includes)
