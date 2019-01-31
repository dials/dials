"""
Tests for the reflection selection algorithm.
"""
import os
import itertools
from cctbx import sgtbx
from dxtbx.serialize import load
from dials.array_family import flex
from dials.algorithms.scaling.Ih_table import IhTable
from dials.algorithms.scaling.reflection_selection import \
  select_highly_connected_reflections,\
  select_connected_reflections_across_datasets,\
  select_highly_connected_reflections_in_bin

def test_select_highly_connected_reflections_in_bin():
  """Test the single-bin selection algorithm."""
  r1 = flex.reflection_table()
  n_list = [3, 3, 2, 1, 1, 2, 2]
  miller_indices = [[(0, 0, i+1)]*n for i, n in enumerate(n_list)]
  r1['miller_index'] = flex.miller_index(
    list(itertools.chain.from_iterable(miller_indices)))
  r1['class_index'] = flex.int([0, 1, 1, 0, 1, 2, 0, 0, 2, 1, 1, 2, 0, 1])
  r1['intensity'] = flex.double(sum(n_list), 1)
  r1['variance'] = flex.double(sum(n_list), 1)
  r1['inverse_scale_factor'] = flex.double(sum(n_list), 1)

  sg = sgtbx.space_group('P1')

  indices, total_in_classes = select_highly_connected_reflections_in_bin(
    r1, sg, min_per_class=2, min_total=6, max_total=100)
  assert list(total_in_classes) == [2, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0]
  assert list(indices) == [0, 1, 2, 3, 4, 5, 10, 11]

def test_select_connected_reflections_across_datasets():
  """Test the basic cross-dataset reflection selection algorithm.

  Make three reflection tables with the following reflections:
             symmetry groups
             0  1  2  3  4  5  6
          0  3  3  2  0  1  1  1
  classes 1  0  2  0  0  3  2  1
          2  2  1  1  5  0  4  0

  With target=5, expect:
  number of chosen reflections per class: [8, 7, 7]
  symmetry groups used:                   [1, 5, 0, 4]
  """

  n1 = [3, 3, 2, 0, 1, 1, 1]
  n2 = [0, 2, 0, 0, 3, 2, 1]
  n3 = [2, 1, 1, 5, 0, 4, 0]

  def make_refl_table(n_list, class_idx=0):
    """Make a reflection table with groups based on n_list."""
    r1 = flex.reflection_table()
    miller_indices = [[(0, 0, i+1)]*n for i, n in enumerate(n_list)]
    print(miller_indices)
    r1['miller_index'] = flex.miller_index(
      list(itertools.chain.from_iterable(miller_indices)))
    r1['intensity'] = flex.double(sum(n_list), 1)
    r1['variance'] = flex.double(sum(n_list), 1)
    r1['inverse_scale_factor'] = flex.double(sum(n_list), 1)
    r1['class_index'] = flex.int(sum(n_list), class_idx)
    return r1

  reflections = [make_refl_table(n1, 0), make_refl_table(n2, 1),
    make_refl_table(n3, 2)]

  space_group = sgtbx.space_group('P1')
  table = IhTable(reflections, space_group)
  indices, datset_ids, total_in_classes = \
    select_connected_reflections_across_datasets(
    table, min_per_class=5, Isigma_cutoff=0.0)
  assert list(total_in_classes) == [8, 7, 7]
  assert list(indices) == [
    0, 1, 2, 3, 4, 5, 8, 9,
    0, 1, 2, 3, 4, 5, 6,
    0, 1, 2, 9, 10, 11, 12]
  assert list(datset_ids) == [0] * 8 + [1] * 7 + [2] * 7

def test_reflection_selection(dials_regression):
  """Use a real dataset to test the selection algorithm."""
  data_dir = os.path.join(dials_regression, "xia2-28",)
  pickle_path = os.path.join(data_dir, "20_integrated.pickle")
  sweep_path = os.path.join(data_dir, "20_integrated_experiments.json")
  reflection_table = flex.reflection_table.from_pickle(pickle_path)
  experiment = load.experiment_list(sweep_path, check_format=False)[0]

  reflection_table['intensity'] = reflection_table['intensity.sum.value']
  reflection_table['variance'] = reflection_table['intensity.sum.variance']
  reflection_table['inverse_scale_factor'] = flex.double(
    reflection_table.size(), 1.0)
  reflection_table = reflection_table.select(reflection_table['variance'] > 0)
  reflection_table = reflection_table.select(reflection_table.get_flags(
    reflection_table.flags.integrated, all=True))

  indices = select_highly_connected_reflections(reflection_table, experiment,
    min_per_area=10, n_resolution_bins=10)
  assert len(indices) > 1650 and len(indices) < 1800
  ##FIXME ^^ different answer on linux/windows? 1658 vs 1712?

  # Give a high min_per_area to check that all reflections with multiplciity > 1
  # are selected.
  indices = select_highly_connected_reflections(reflection_table, experiment,
    min_per_area=50, n_resolution_bins=10)
  # this dataset has 48 reflections with multiplicity = 1
  assert len(indices) == reflection_table.size() - 48
