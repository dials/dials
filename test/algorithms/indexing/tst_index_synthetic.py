from __future__ import division

import math
import os
import random
import libtbx.load_env
from libtbx import easy_pickle
from libtbx.test_utils import open_tmp_directory
from scitbx import matrix
from cctbx import crystal, sgtbx, uctbx

from dxtbx.model.experiment.experiment_list import Experiment
from dxtbx.serialize import dump, load

from dials.array_family import flex


have_dials_regression = libtbx.env.has_module("dials_regression")
if have_dials_regression:
  dials_regression = libtbx.env.find_in_repositories(
    relative_path="dials_regression",
    test=os.path.isdir)


def random_rotation():
  import random
  from scitbx.math import euler_angles_as_matrix
  return euler_angles_as_matrix([random.uniform(0,360) for i in xrange(3)])


def any_compatible_unit_cell(space_group, volume=None, asu_volume=None):
  # based on sgtbx.any_compatible_unit_cell()
  assert [volume, asu_volume].count(None) == 1
  if volume is None:
    volume = asu_volume * space_group.order_z()
  sg_number = space_group.type().number()
  a = 1.
  b = random.uniform(0.4, 3)
  c = random.uniform(0.4, 3)
  while True:
    alpha = random.uniform(0, 180)
    beta = random.uniform(0, 180)
    gamma = random.uniform(0, 180)

    # Allowed combinations of angles for a triclinic cell:
    #   http://dx.doi.org/10.1107/S0108767310044296
    if ((alpha + beta + gamma) > 0 and
        (alpha + beta + gamma) < 360 and
        (alpha + beta - gamma) > 0 and
        (alpha + beta - gamma) < 360 and
        (alpha - beta + gamma) > 0 and
        (alpha - beta + gamma) < 360 and
        (-alpha + beta + gamma) > 0 and
        (-alpha + beta + gamma) < 360):
      break
    else:
      print a, b, c, alpha, beta, gamma

  if sg_number <   3:
    params = (a, b, c, alpha, beta, gamma)
  elif sg_number <  16:
    params = (a, b, c, 90, beta, 90)
  elif sg_number <  75:
    params = (a, b, c, 90, 90, 90)
  elif sg_number < 143:
    params = (a, a, c, 90, 90, 90)
  elif sg_number < 195:
    params = (a, a, c, 90, 90, 120)
  else:
    params = (a, a, a, 90, 90, 90)
  unit_cell = uctbx.unit_cell(params).change_basis(
    cb_op=space_group.info().change_of_basis_op_to_reference_setting().inverse())
  f = (volume / unit_cell.volume())**(1/3.)
  params = list(unit_cell.parameters())
  for i in xrange(3): params[i] *= f
  return uctbx.unit_cell(params)


def generate_spots(crystal_model, detector, beam, goniometer=None, scan=None,
                   sel_fraction=1.0):
  import math

  experiment = Experiment(beam=beam,
                          detector=detector,
                          goniometer=goniometer,
                          scan=scan,
                          crystal=crystal_model)

  # if we don't set the imageset then from_predictions uses the StillsReflectionPredictor :-(
  from dxtbx.imageset import NullReader, ImageSweep
  imageset = ImageSweep(NullReader, indices=range(len(scan.get_epochs())), beam=beam, goniometer=goniometer,
                        detector=detector, scan=scan)
  experiment.imageset = imageset

  predicted = flex.reflection_table.from_predictions(experiment)

  sel = flex.random_selection(len(predicted),
                              int(math.floor(sel_fraction*len(predicted))))
  predicted = predicted.select(sel)
  predicted['imageset_id'] = flex.size_t(len(predicted), 0)
  predicted['xyzobs.px.value'] = predicted['xyzcal.px']
  predicted['xyzobs.px.variance'] = flex.vec3_double(
    len(predicted), (0.5,0.5,0.5))
  return predicted


def generate_crystal(unit_cell, space_group):
  from dxtbx.model.crystal import crystal_model

  B = matrix.sqr(unit_cell.fractionalization_matrix()).transpose()
  U = random_rotation()
  A = U * B
  direct_matrix = A.inverse()
  return crystal_model(direct_matrix[0:3],
                       direct_matrix[3:6],
                       direct_matrix[6:9],
                       space_group=space_group)


def run_indexing(datablock, strong_spots, crystal_model, rmsds):

  cwd = os.path.abspath(os.curdir)
  tmp_dir = os.path.abspath(open_tmp_directory(suffix="test_dials_index"))
  os.chdir(tmp_dir)

  sweep_path = os.path.join(tmp_dir, "datablock.json")
  pickle_path = os.path.join(tmp_dir, "strong.pickle")

  dump.datablock(datablock, sweep_path)
  easy_pickle.dump(pickle_path, strong_spots)

  from dials.test.algorithms.indexing.tst_index import run_one_indexing

  space_group_info = crystal_model.get_space_group()
  symmetry = crystal.symmetry(unit_cell=crystal_model.get_unit_cell(),
                              space_group=crystal_model.get_space_group())

  expected_rmsds = [1.1*r for r in rmsds]

  imageset = datablock[0].extract_imagesets()[0]
  pixel_size = imageset.get_detector()[0].get_pixel_size()
  phi_width = imageset.get_scan().get_oscillation()[1] * math.pi/180

  expected_rmsds = [1.1 * rmsds[0] * pixel_size[0],
                    1.1 * rmsds[1] * pixel_size[1],
                    1.1 * rmsds[2] * phi_width]

  run_one_indexing(pickle_path=pickle_path, sweep_path=sweep_path,
                   extra_args=[],
                   expected_unit_cell=symmetry.minimum_cell().unit_cell(),
                   expected_rmsds=expected_rmsds,
                   #expected_hall_symbol=crystal_model.get_space_group().type().hall_symbol(),
                   expected_hall_symbol=' P 1',
                   )


def add_random_noise_xyz(datablock, strong_spots, rmsds):
  import scitbx.random
  errors = flex.vec3_double(*list(
    scitbx.random.variate(
      scitbx.random.normal_distribution(0, rmsds[i]))(len(strong_spots))
    for i in range(3)))
  strong_spots['xyzobs.px.value'] += errors


def run():

  #space_group = sgtbx.space_group_info(number=random.randint(1,230)).group()
  space_group = sgtbx.space_group_info(
    number=random.choice((1,3,16,75,143,195))).group()

  unit_cell = any_compatible_unit_cell(
    space_group, volume=random.uniform(1e3, 1e6))

  fname = os.path.join(dials_regression, "centroid_test_data", "datablock.json")
  datablock = load.datablock(fname, check_format=False)
  imageset = datablock[0].extract_imagesets()[0]
  scan = imageset.get_scan()
  scan.set_image_range((1,900))

  crystal_model = generate_crystal(unit_cell, space_group)
  print crystal_model
  print unit_cell.minimum_cell()
  strong_spots = generate_spots(
    crystal_model, imageset.get_detector(),
    imageset.get_beam(),
    goniometer=imageset.get_goniometer(),
    scan=scan,
    sel_fraction=0.25)

  rmsds = (0.5, 0.5, 0.5) # px/image
  #print strong_spots[0]['xyzobs.px.value']
  add_random_noise_xyz(datablock, strong_spots, rmsds)
  #print strong_spots[0]['xyzobs.px.value']

  run_indexing(datablock, strong_spots, crystal_model, rmsds)

  print "OK"

if __name__ == '__main__':
  run()
