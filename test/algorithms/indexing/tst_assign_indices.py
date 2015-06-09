from __future__ import division

import os
import libtbx.load_env
from cctbx import crystal, miller, sgtbx
from scitbx import matrix
from dxtbx.serialize import load
from dxtbx.model.experiment.experiment_list import Experiment, ExperimentList
from dxtbx.model.crystal import crystal_model
from dials.array_family import flex

# set random seeds so tests more reliable
seed = 54321
import random
random.seed(seed)
flex.set_random_seed(seed)

have_dials_regression = libtbx.env.has_module("dials_regression")
if have_dials_regression:
  dials_regression = libtbx.env.find_in_repositories(
    relative_path="dials_regression",
    test=os.path.isdir)
else:
  print 'SKIP: dials_regression not configured'
  exit(0)


def random_rotation(angle_min=0, angle_max=360):
  import random
  from scitbx.math import euler_angles_as_matrix
  return euler_angles_as_matrix(
    [random.uniform(angle_min,angle_max) for i in xrange(3)])


def run(space_group_info):
  datablock_json = os.path.join(
    dials_regression, "indexing_test_data",
    "i04_weak_data", "datablock_orig.json")

  datablock = load.datablock(datablock_json, check_format=False)[0]
  sweep = datablock.extract_imagesets()[0]

  sweep._indices = sweep._indices[:20]
  sweep.set_scan(sweep.get_scan()[:20])

  import random
  space_group = space_group_info.group()
  unit_cell = space_group_info.any_compatible_unit_cell(volume=random.uniform(1e4,1e6))

  crystal_symmetry = crystal.symmetry(unit_cell=unit_cell,
                                      space_group=space_group)
  crystal_symmetry.show_summary()

  # the reciprocal matrix
  B = matrix.sqr(unit_cell.fractionalization_matrix()).transpose()
  U = random_rotation()
  A = U * B

  direct_matrix = A.inverse()
  cryst_model = crystal_model(direct_matrix[0:3],
                              direct_matrix[3:6],
                              direct_matrix[6:9],
                              space_group=space_group)
  experiment = Experiment(imageset=sweep,
                          beam=sweep.get_beam(),
                          detector=sweep.get_detector(),
                          goniometer=sweep.get_goniometer(),
                          scan=sweep.get_scan(),
                          crystal=cryst_model)
  predicted_reflections = flex.reflection_table.from_predictions(
    experiment)
  use_fraction = 0.3
  use_sel = flex.random_selection(
    len(predicted_reflections), int(use_fraction*len(predicted_reflections)))
  predicted_reflections = predicted_reflections.select(use_sel)
  miller_indices = predicted_reflections['miller_index']
  miller_set = miller.set(
    crystal_symmetry, miller_indices, anomalous_flag=True)
  predicted_reflections['xyzobs.mm.value'] = predicted_reflections['xyzcal.mm']
  predicted_reflections['id'] = flex.size_t(len(predicted_reflections), 0)
  from dials.algorithms.indexing.indexer import indexer_base
  indexer_base.map_centroids_to_reciprocal_space(
    predicted_reflections, sweep.get_detector(), sweep.get_beam(),
    sweep.get_goniometer())


  # check that local and global indexing worked equally well in absence of errors
  result = compare_global_local(experiment, predicted_reflections,
                                miller_indices)
  assert result.misindexed_local == 0
  assert result.misindexed_global == 0

  a, b, c = cryst_model.get_real_space_vectors()
  relative_error = 0.02
  a *= (1+relative_error)
  b *= (1+relative_error)
  c *= (1+relative_error)

  cryst_model2 = crystal_model(a, b, c, space_group=space_group)
  experiment.crystal = cryst_model2

  result = compare_global_local(experiment, predicted_reflections,
                                miller_indices)

  # check that the local indexing did a better job given the errors in the basis vectors
  #assert result.misindexed_local < result.misindexed_global
  assert result.misindexed_local == 0
  assert result.correct_local > result.correct_global
  # usually the number misindexed is much smaller than this
  assert result.misindexed_local < (0.001 * len(result.reflections_local))

  # the reciprocal matrix
  A = cryst_model.get_A()
  A = random_rotation(angle_max=0.03) * A

  direct_matrix = A.inverse()
  cryst_model2 = crystal_model(direct_matrix[0:3],
                               direct_matrix[3:6],
                               direct_matrix[6:9],
                               space_group=space_group)
  experiment.crystal = cryst_model2

  result = compare_global_local(experiment, predicted_reflections,
                                miller_indices)

  # check that the local indexing did a better job given the errors in the basis vectors
  assert result.misindexed_local <= result.misindexed_global, (
    result.misindexed_local, result.misindexed_global)
  assert result.misindexed_local < 0.01 * result.correct_local
  assert result.correct_local > result.correct_global
  # usually the number misindexed is much smaller than this
  assert result.misindexed_local < (0.001 * len(result.reflections_local))





class compare_global_local(object):

  def __init__(self, experiment, reflections,
               expected_miller_indices):

    from dials.algorithms.indexing \
         import index_reflections, index_reflections_local
    import copy

    # index reflections using simple "global" method
    self.reflections_global = copy.deepcopy(reflections)
    self.reflections_global['id'] = flex.int(len(self.reflections_global), -1)
    self.reflections_global['imageset_id'] = flex.int(len(self.reflections_global), 0)
    index_reflections(
      self.reflections_global, ExperimentList([experiment]))
    non_zero_sel = (self.reflections_global['miller_index'] != (0,0,0))
    assert self.reflections_global['id'].select(~non_zero_sel).all_eq(-1)
    self.misindexed_global = (
      expected_miller_indices == self.reflections_global['miller_index']).select(
        non_zero_sel).count(False)
    self.correct_global = (
      expected_miller_indices == self.reflections_global['miller_index']).count(True)


    # index reflections using xds-style "local" method
    self.reflections_local = copy.deepcopy(reflections)
    self.reflections_local['id'] = flex.int(len(self.reflections_local), -1)
    index_reflections_local(
      self.reflections_local, ExperimentList([experiment]))
    non_zero_sel = (self.reflections_local['miller_index'] != (0,0,0))
    assert self.reflections_local['id'].select(~non_zero_sel).all_eq(-1)
    self.misindexed_local = (
      expected_miller_indices == self.reflections_local['miller_index']).select(
        non_zero_sel).count(False)
    self.correct_local = (
      expected_miller_indices == self.reflections_local['miller_index']).count(True)

    print self.misindexed_global, self.correct_global, len(self.reflections_global)
    print self.misindexed_local, self.correct_local, len(self.reflections_local)


if __name__ == '__main__':
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  from cctbx.sgtbx import bravais_types
  for symbol in bravais_types.acentric:
    space_group_info = sgtbx.space_group_info(symbol=symbol)
    run(space_group_info)
