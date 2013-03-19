from __future__ import division
def fractional_difference(target, reference):
    import math
    return math.fabs(target - reference) / reference

def angular_difference(target, reference):
    import math

    angle = target.angle(reference)

    if angle < math.pi - angle:
        return angle
    else:
        return math.pi - angle

def tst_toy_centroid_helpers():
    from dials.algorithms.centroid.toy_centroid_helpers import \
        generate_fake_profile, form_covariance_matrix
    from scitbx import matrix

    x = matrix.col((1, 1, 0)).normalize()
    y = matrix.col((1, -1, 0)).normalize()
    z = matrix.col((0, 0, 1))

    xyz = (x, y, z)
    variances = (4, 2, 1)

    pixels = generate_fake_profile(xyz, variances, 100000)

    evalues, evectors = form_covariance_matrix(pixels, 0, 0, 0)

    for j in range(3):
        assert(fractional_difference(evalues[j], variances[j]) < 0.1)
        assert(angular_difference(evectors[j], xyz[j]) < 0.1)

    print 'OK'

if __name__ == '__main__':
    tst_toy_centroid_helpers()
