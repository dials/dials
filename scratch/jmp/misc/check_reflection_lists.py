

def equal(a, b):
    try:
        assert(a == b)
    except Exception, e:
        print a, b
        raise e

def almost_equal(a, b, eps=1e-7):
    try:
        assert(abs(a - b) < eps)
    except Exception, e:
        print a, b
        raise e

import sys
from dials.model.data import ReflectionList

filename1 = sys.argv[1]
filename2 = sys.argv[2]

from dials.model.serialize import load
from cctbx.array_family import flex

refl1 = load.reflections(filename1)
refl2 = load.reflections(filename2)

print "Length: ", len(refl1), len(refl2)
#assert(len(refl1) == len(refl2))

refl1 = ReflectionList([r for r in refl1 if r.is_valid()])
refl2 = ReflectionList([r for r in refl2 if r.is_valid()])

num_valid1 = len([r for r in refl1 if r.is_valid()])
num_valid2 = len([r for r in refl2 if r.is_valid()])
print "Valid: ", num_valid1, num_valid2
assert(num_valid1 == num_valid2)
from scitbx import matrix
print "Checking"
for r1, r2 in zip(refl1, refl2):
    equal(r1.is_valid(), r2.is_valid())
    equal(r1.miller_index, r2.miller_index)
    almost_equal(r1.rotation_angle, r2.rotation_angle)
    almost_equal(matrix.col(r1.beam_vector), matrix.col(r2.beam_vector))
    almost_equal(matrix.col(r1.image_coord_px), matrix.col(r2.image_coord_px))
    almost_equal(r1.frame_number, r2.frame_number)
    almost_equal(matrix.col(r1.centroid_position), matrix.col(r2.centroid_position))
    almost_equal(r1.intensity, r2.intensity)
    almost_equal(r1.intensity_variance, r2.intensity_variance)
    almost_equal(r1.corrected_intensity, r2.corrected_intensity)
    almost_equal(r1.corrected_intensity_variance, r2.corrected_intensity_variance)
