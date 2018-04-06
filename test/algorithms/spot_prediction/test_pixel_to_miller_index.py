from __future__ import absolute_import, division, print_function

import os
import pytest

def test(dials_regression):
  from dials.algorithms.spot_prediction import PixelToMillerIndex
  from dials.array_family import flex
  from dxtbx.model.experiment_list import ExperimentListFactory

  filename = os.path.join(dials_regression, "centroid_test_data", "experiments.json")

  experiments = ExperimentListFactory.from_json_file(filename)

  transform = PixelToMillerIndex(
    experiments[0].beam,
    experiments[0].detector,
    experiments[0].goniometer,
    experiments[0].scan,
    experiments[0].crystal)

  reflections = flex.reflection_table.from_predictions(experiments[0])

  for r in reflections:
    panel = r['panel']
    x, y, z = r['xyzcal.px']
    h0 = r['miller_index']
    h1 = transform.h(panel, x, y, z)
    assert h0 == pytest.approx(h1, abs=1e-7)
