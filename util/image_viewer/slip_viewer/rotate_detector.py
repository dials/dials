from __future__ import division, print_function
import sys
from dxtbx.model.experiment_list import ExperimentListFactory
from scitbx.matrix import col

"""
Example script for rotating a detector in an experiment.json
"""

experiments = ExperimentListFactory.from_json_file(sys.argv[1])
print("Original detector")
detector = experiments[0].detector
print(detector)

h = detector.hierarchy()

fast = col(h.get_fast_axis())
slow = col(h.get_slow_axis())
origin = col(h.get_origin())

normal = fast.cross(slow)
rotation = normal.axis_and_angle_as_r3_rotation_matrix(90, deg=True)

h.set_frame(rotation * fast, rotation * slow, rotation * origin)

print("Rotated detector")
print(detector)

experiments.as_json("rotated.expt")
