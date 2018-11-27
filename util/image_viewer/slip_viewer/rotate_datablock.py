from __future__ import division, print_function
import sys
from dxtbx.datablock import DataBlockFactory, DataBlockDumper
from scitbx.matrix import col

"""
Example script for rotating a detector in a datablock.json
"""

datablock = DataBlockFactory.from_json_file(sys.argv[1])[0]
print("Original detector")
detector = datablock.unique_detectors()[0]
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

dump = DataBlockDumper(datablock)
dump.as_json("rotated.json")
