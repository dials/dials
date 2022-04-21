from __future__ import annotations

import sys

from dxtbx.model.experiment_list import ExperimentList
from scitbx.matrix import col

"""
Example script for rotating a detector in an experiment.json
"""

experiments = ExperimentList.from_file(sys.argv[1])
print("Original detector")
detector = experiments[0].detector
print(detector)

h = detector.hierarchy()

fast = col(h.get_fast_axis())
slow = col(h.get_slow_axis())
origin = col(h.get_origin())

normal = col(h.get_normal())
rotation = normal.axis_and_angle_as_r3_rotation_matrix(90, deg=True)

h.set_frame(rotation * fast, rotation * slow, rotation * origin)

print("Rotated detector")
print(detector)

experiments.as_json("rotated.expt")
