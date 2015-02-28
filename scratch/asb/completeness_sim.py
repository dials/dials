from __future__ import division
#from iotbx import reflection_file_reader
from cctbx.crystal import symmetry
from cctbx.array_family import flex
from scitbx.matrix import sqr,col
from dxtbx.model.detector import detector_factory
from dials.array_family import flex
from cctbx import crystal
from dxtbx.model.crystal import crystal_model
from dxtbx.model.beam import beam_factory
from dxtbx.model.experiment.experiment_list import Experiment
from cctbx import miller
from math import sqrt

# Set up dxtbx models, starting with a detector with a single panel, similar to a CSPAD
#detector=detector_factory.simple('SENSOR_UNKNOWN',100,(97.075,97.075),'+x','-y',(0.11,0.11),(1765,1765))

# Here's the MarCCD at XPP
#detector=detector_factory.simple('SENSOR_UNKNOWN',150,(162.5,162.5),'+x','-y',(0.079346,0.079346),(4096,4096))

# Here we simulate the gaps in a CSPAD by using a CSPAD CBF
import dxtbx, libtbx.load_env, os
img_path = os.path.join(libtbx.env.find_in_repositories("dials_regression"), "spotfinding_test_data/idx-s00-20131106040304531.cbf")
img = dxtbx.load(img_path)
detector = img.get_detector()

# Manually push the detector in to 100 mm
h = detector.hierarchy()
h.set_local_frame(h.get_fast_axis(), h.get_slow_axis(), (0,0,-100))

# Beam model
wavelength=1.32 # from Boutet 2012
beam = beam_factory.simple_directional((0,0,1), wavelength)

# Unit cell and space group for lysozyme from Boutet 2012
known_symmetry=crystal.symmetry("79,79,38,90,90,90","P43212")
#known_symmetry=crystal.symmetry("79,79,38,90,90,90","P1")
uc=known_symmetry.unit_cell()
spgrp=known_symmetry.space_group()

def get_random_predictions():
  """ Return a DIALS reflection table representing predictions using the given models.
  Assumes a Ewald proximity model for mosaicity """
  # The U matrix to calculate A*
  rot = flex.random_double_r3_rotation_matrix()
  A = sqr(rot) * sqr(uc.reciprocal().orthogonalization_matrix())
  A_inv = A.inverse()

  # Use the matrix to create a crystal model
  a = col(A_inv[:3])
  b = col(A_inv[3:6])
  c = col(A_inv[6:])
  cm = crystal_model(a, b, c, space_group=spgrp)

  # This DIALS object has no gonio or scan which will identify it as a still
  expt = Experiment(beam=beam,
                    detector=detector,
                    goniometer=None,
                    scan=None,
                    crystal=cm)

  # Predict the reflections
  return flex.reflection_table.from_predictions(expt)

number_xtals = 0
while True:
  # Add crystals to the miller set until a certain threshold is reached
  number_xtals += 1

  predictions = get_random_predictions()

  # Use cctbx's miller set to handle binning
  if 'miller_set' not in locals():
    miller_set = miller.set(known_symmetry, predictions['miller_index'], anomalous_flag=True)
  else:
    miller_set.indices().extend(predictions['miller_index'])
  miller_set = miller_set.resolution_filter(d_min=2)

  # Report multiplicity
  binner = miller_set.setup_binner(n_bins=10)
  data = []
  for i in binner.range_all():
    if binner.counts_complete()[i] > 0:
      data.append(binner.counts()[i] / binner.counts_complete()[i])
    else:
      data.append(0)

  # Report completeness
  unique_set = miller_set.unique_under_symmetry()
  unique_set.setup_binner(n_bins=10)
  binned_data = unique_set.completeness(use_binning=True,multiplier=100)

  print "Crystals:", number_xtals
  print "Multiplicity"
  binner.show_data(data, data_fmt="%.3f")
  print "Completeness"
  binned_data.show()
  n_obs = len(miller_set.indices())
  print "N obs:", n_obs, "N unique:", len(unique_set.indices())
  n_sites = 1; f = 2; ccano = 0.36
  print "Estimated Sanom assuming CCano %.2f, %d site and f=%.2f: %.3f"%(ccano, n_sites, f, 0.36*sqrt(n_obs)/(sqrt(1)*sqrt(2)))

  # Thresholding
  if [binned_data.data[i] < 95 for i in binned_data.binner.range_used()].count(False) != 10:
    continue

  if [binner.counts()[i] / binner.counts_complete()[i] < 4 for i in binner.range_used()].count(False) == 10:
    break

print "Took", number_xtals, "crystals to 95% completeness, 4 fold redundancy"
