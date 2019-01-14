"""
Module of functions for handling operations on lists of reflection tables
and experiment lists.
"""

from __future__ import print_function

import logging
from dials.array_family import flex
from libtbx.utils import Sorry

logger = logging.getLogger('dials')

def parse_multiple_datasets(reflections):
  """
  Split a list of multi-dataset reflection tables, selecting on id

  If duplicate id values are found, the id columns are renumbered from 0..n-1,
  taking care of experiment identifiers if these are set.

  Args:
      reflections (list): a list of reflection tables, each of which may contain
          multiple datasets

  Returns:
      (list): a list of reflection tables corresponding to single datasets

  """
  single_reflection_tables = []
  dataset_id_list = []
  for refl_table in reflections:
    dataset_ids = set(refl_table['id']).difference(set([-1]))
    dataset_id_list.extend(list(dataset_ids))
    if len(dataset_ids) > 1:
      logger.info('Detected existence of a multi-dataset reflection table \n'
        'containing %s datasets. \n', len(dataset_ids))
      ##FIXME fix split_by_experiment_id so that don't need to filter
      # unindxeded reflections here to get rid of id = -1
      if -1 in refl_table['id']:
        refl_table = refl_table.select(refl_table['id'] != -1)
      result = refl_table.split_by_experiment_id()
      single_reflection_tables.extend(result)
    else:
      single_reflection_tables.append(refl_table)
  if len(dataset_id_list) != len(set(dataset_id_list)): # need to reset some ids
    logger.warn('Duplicate dataset ids found in different reflection tables. \n'
    'These will be treated as coming from separate datasets, and \n'
    'new dataset ids will be assigned for the whole dataset. \n')
    new_id_list = range(len(dataset_id_list))
    for r, old_id, new_id in zip(single_reflection_tables, dataset_id_list,
      new_id_list):
      r['id'] = flex.int(r.size(), new_id)
      if list(r.experiment_identifiers()):# if identifiers, need to update
        expid = r.experiment_identifiers()[old_id]
        del r.experiment_identifiers()[old_id]
        r.experiment_identifiers()[new_id] = expid
  return single_reflection_tables

def get_next_unique_id(unique_id, used_ids):
  """
  Test a list of used id strings to see if it contains str(unique_id)

  Args:
      unique_id (int): Integer value to be converted to str
      used_ids (list): A list of strings

  Returns:
      (int): The lowest int >= unique_id for which str(int) is not in used_ids

  """
  while str(unique_id) in used_ids:
    unique_id += 1
  return unique_id

def assign_unique_identifiers(experiments, reflections, identifiers=None):
  """
  Assign unique experiment identifiers to experiments and reflections lists.

  If experiment identifiers are not set for some datasets, then new unique
  identifiers are given to those, and the 'id' column for all reflection tables
  are set sequentially from 0..n-1.

  Args:
      experiments: An ExperimentList
      reflections (list): A list of reflection tables

  Returns:
      (tuple): tuple containing:
          experiments: The updated ExperimentList
          reflections (list): A list of the updated reflection tables

  Raises:
      Sorry: If the number of reflection tables and experiments are unequal,
          after attempted parsing of multi-dataset reflection tables.
  """
  if len(experiments) != len(reflections):
    reflections = parse_multiple_datasets(reflections)
    if len(reflections) != len(experiments):
      raise Sorry("""Unable to split a list of reflection tables to be the same
length as the experiments list. Please check the input data.""")
  if identifiers:
    if len(identifiers) != len(reflections):
      raise Sorry("""Number of provided identifiers (%s) not the same length as
number of datasets (%s)""" % (len(identifiers), len(reflections)))
    for i, (exp, refl) in enumerate(zip(experiments, reflections)):
      exp.identifier = identifiers[i]
      for k in refl.experiment_identifiers().keys():
        del refl.experiment_identifiers()[k]
      refl.experiment_identifiers()[i] = identifiers[i]
      refl['id'] = flex.int(refl.size(), i)
  used_str_ids = []
  for exp, refl in zip(experiments, reflections):
    if exp.identifier != '':
      assert list(refl.experiment_identifiers().values()) == [exp.identifier]
      used_str_ids.append(exp.identifier)
  if len(set(used_str_ids)) == len(reflections): #all set, don't do anything
    pass
  else: #keep identifiers if set, reset table id column from 0..n-1
    unique_id = 0
    for i, (exp, refl) in enumerate(zip(experiments, reflections)):
      if exp.identifier == '':
        unique_id = get_next_unique_id(unique_id, used_str_ids)
        strid = '%i' % unique_id
        exp.identifier = strid
        refl.experiment_identifiers()[i] = strid
        unique_id += 1
      else:
        k = refl.experiment_identifiers().keys()[0]
        expid = refl.experiment_identifiers().values()[0]
        del refl.experiment_identifiers()[k]
        refl.experiment_identifiers()[i] = expid
      refl['id'] = flex.int(refl.size(), i)
  return experiments, reflections

def select_datasets_on_ids(experiments, reflection_table_list,
  exclude_datasets=None, use_datasets=None):
  """
  Select a subset of experiments and reflection tables based on identifiers.

  This performs a similar function to the select/remove_on_experiment_identifiers
  methods of ExperimentList and reflection_table, with additional logic to handle
  the case of a list of reflection tables, rather than a single one.

  Args:
      experiments: An ExperimentList
      reflection_table_list (list): a list of reflection tables
      exclude_datasets (list): a list of experiment_identifiers to exclude
      use_datasets (list): a list of experiment_identifiers to use

  Returns:
      (tuple): tuple containing:
          experiments: The updated ExperimentList
          list_of_reflections (list): A list of the updated reflection tables

  Raises:
      Sorry: If both use_datasets and exclude datasets are used, if not all
          experiment identifiers are set, if an identifier in exclude_datasets
          or use_datasets is not in the list.
  """
  if not use_datasets and not exclude_datasets:
    return experiments, reflection_table_list
  if use_datasets and exclude_datasets:
    raise Sorry("The options use_datasets and exclude_datasets cannot be used in conjuction.")
  if experiments.identifiers().count('') > 0:
    logger.warn('\nERROR: Attempting to choose datasets based on unique identifier,\n'
      'but not all datasets currently have a unique identifier! Please make\n'
      'sure all identifiers are set before attempting to select datasets.\n')
    logger.info('Current identifiers set as: %s', list(experiments.identifiers()))
    raise Sorry("Not all experiment identifiers set in the ExperimentList")
  list_of_reflections = []
  if use_datasets:
    total_found = 0
    for reflection_table in reflection_table_list:
      expids = list(reflection_table.experiment_identifiers().values())
      expids_in_this_table = list(set(use_datasets).intersection(set(expids)))
      total_found += len(expids_in_this_table)
      if expids_in_this_table: #if none, then no datasets wanted from table
        list_of_reflections.append(reflection_table.select_on_experiment_identifiers(
          expids_in_this_table))
    if total_found != len(use_datasets):
      raise Sorry("""Attempting to select datasets based on identifiers that
are not found in the experiment list / reflection tables.""")
    experiments.select_on_experiment_identifiers(use_datasets)
  elif exclude_datasets:
    if len(set(exclude_datasets).intersection(set(experiments.identifiers()))) != \
      len(set(exclude_datasets)):
      raise Sorry("""Attempting to exclude datasets based on identifiers that
are not found in the experiment list / reflection tables.""")
    for reflection_table in reflection_table_list:
      expids = list(reflection_table.experiment_identifiers().values())
      expids_in_this_table = list(set(exclude_datasets).intersection(set(expids)))
      if expids_in_this_table: #only append if data left after removing
        r_table = reflection_table.remove_on_experiment_identifiers(
          expids_in_this_table)
        if r_table.size() > 0:
          list_of_reflections.append(r_table)
      else:
        list_of_reflections.append(reflection_table)
    experiments.remove_on_experiment_identifiers(exclude_datasets)
  return experiments, list_of_reflections
