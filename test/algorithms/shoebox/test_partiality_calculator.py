from __future__ import absolute_import, division, print_function

import math
import os

def test(dials_regression):
  filename = os.path.join(dials_regression, 'centroid_test_data', 'fake_long_experiments.json')

  from dxtbx.model.experiment_list import ExperimentListFactory
  from dxtbx.model.experiment_list import ExperimentList
  exlist = ExperimentListFactory.from_json_file(filename)
  assert len(exlist) == 1
  experiment = exlist[0]

  # Set the delta_divergence/mosaicity
  n_sigma = 5
  sigma_b = 0.060 * math.pi / 180
  sigma_m = 0.154 * math.pi / 180

  from dials.algorithms.profile_model.gaussian_rs import Model
  profile_model = Model(None, n_sigma, sigma_b, sigma_m)
  experiment.profile = profile_model
  experiments = ExperimentList()
  experiments.append(experiment)

  from dials.algorithms.profile_model.gaussian_rs import \
    PartialityCalculator3D
  from dials.array_family import flex

  calculator = PartialityCalculator3D(
    experiment.beam,
    experiment.goniometer,
    experiment.scan,
    sigma_m)

  predicted = flex.reflection_table.from_predictions_multi(experiments)
  predicted['bbox'] = predicted.compute_bbox(experiments)

  # Remove any touching edges of scan to get only fully recorded
  x0, x1, y0, y1, z0, z1 = predicted['bbox'].parts()
  predicted = predicted.select((z0 > 0) & (z1 < 100))
  assert(len(predicted) > 0)

  # Compute partiality
  partiality = calculator(
    predicted['s1'],
    predicted['xyzcal.px'].parts()[2],
    predicted['bbox'])

  # Should have all fully recorded
  assert len(partiality) == len(predicted)
  three_sigma = 0.5 * (math.erf(3.0 / math.sqrt(2.0)) - math.erf(-3.0 / math.sqrt(2.0)))
  assert partiality.all_gt(three_sigma)

  # Trim bounding boxes
  x0, x1, y0, y1, z0, z1 = predicted['bbox'].parts()
  z0 = z0 + 1
  z1 = z1 - 1
  predicted['bbox'] = flex.int6(x0, x1, y0, y1, z0, z1)
  predicted = predicted.select(z1 > z0)
  assert len(predicted) > 0

  # Compute partiality
  partiality = calculator(
    predicted['s1'],
    predicted['xyzcal.px'].parts()[2],
    predicted['bbox'])

  # Should have all partials
  assert len(partiality) == len(predicted)
  assert partiality.all_le(1.0) and partiality.all_gt(0)
