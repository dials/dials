from __future__ import division

import math
from scitbx import matrix
from dials.scratch.dgw.prediction import angle_predictor, angle_predictor_py
from cctbx.array_family import flex


#### For a unit test, we need two of Graeme's functions taken from
#### use_case_xds_method/tdi.py
def orthogonal_component(reference, changing):
    '''Return unit vector corresponding to component of changing orthogonal to
    reference.'''

    r = reference.normalize()
    c = changing.normalize()

    return (c - c.dot(r) * r).normalize()

def align_reference_frame(primary_axis, primary_target,
                          secondary_axis, secondary_target):
    '''Compute a rotation matrix R: R x primary_axis = primary_target and
    R x secondary_axis places the secondary_axis in the plane perpendicular
    to the primary_target, as close as possible to the secondary_target.
    Require: primary_target orthogonal to secondary_target, primary axis
    not colinear with secondary axis.'''

    from scitbx import matrix
    import math

    if type(primary_axis) == type(()) or type(primary_axis) == type([]):
        primary_axis = matrix.col(primary_axis).normalize()
    else:
        primary_axis = primary_axis.normalize()

    if type(primary_target) == type(()) or type(primary_target) == type([]):
        primary_target = matrix.col(primary_target).normalize()
    else:
        primary_target = primary_target.normalize()

    if type(secondary_axis) == type(()) or type(secondary_axis) == type([]):
        secondary_axis = matrix.col(secondary_axis).normalize()
    else:
        secondary_axis = secondary_axis.normalize()

    if type(secondary_target) == type(()) or \
           type(secondary_target) == type([]):
        secondary_target = matrix.col(secondary_target).normalize()
    else:
        secondary_target = secondary_target.normalize()

    # check properties of input axes

    assert(math.fabs(primary_axis.angle(secondary_axis) % math.pi) > 0.001)
    assert(primary_target.dot(secondary_target) < 0.001)

    if primary_target.angle(primary_axis) % math.pi:
        axis_p = primary_target.cross(primary_axis)
        angle_p = - primary_target.angle(primary_axis)
        Rprimary = axis_p.axis_and_angle_as_r3_rotation_matrix(angle_p)
    elif primary_target.angle(primary_axis) < 0:
        axis_p = primary_axis.ortho().normalize()
        angle_p = math.pi
        Rprimary = axis_p.axis_and_angle_as_r3_rotation_matrix(angle_p)
    else:
        Rprimary = matrix.identity(3)

    axis_r = secondary_target.cross(Rprimary * secondary_axis)
    axis_s = primary_target
    if (axis_r.angle(primary_target) > 0.5 * math.pi):
        angle_s = orthogonal_component(axis_s, secondary_target).angle(
            orthogonal_component(axis_s, Rprimary * secondary_axis))
    else:
        angle_s = - orthogonal_component(axis_s, secondary_target).angle(
            orthogonal_component(axis_s, Rprimary * secondary_axis))

    Rsecondary = axis_s.axis_and_angle_as_r3_rotation_matrix(angle_s)

    return Rsecondary * Rprimary

if __name__ == "__main__":

    ### Perform a unit test by comparing the results of angle_predictor_py to 
    ### those from rstbx's rotation_angles

    ### python and cctbx imports
    import random
    from rstbx.diffraction import rotation_angles
    from rstbx.diffraction import full_sphere_indices
    from cctbx.sgtbx import space_group, space_group_symbols
    from libtbx.test_utils import approx_equal

    ### import models
    from dials.scratch.dgw.crystal_model import crystal
    from dials.model.experiment import beam_factory, goniometer_factory

    ### local functions
    def random_direction_close_to(vector):
        return vector.rotate_around_origin(matrix.col(
                    (random.random(),
                     random.random(),
                     random.random())).normalize(),
                     random.gauss(0, 1.0),  deg = True)

    ### Create models

    # make a beam vector close to the Z axis
    direction = random_direction_close_to(matrix.col((0, 0, 1)))
    mybeam = beam_factory.make_beam(direction, 1.5)

    # make a random crystal
    a = random.uniform(10,50) * random_direction_close_to(matrix.col((1, 0, 0)))
    b = random.uniform(10,50) * random_direction_close_to(matrix.col((0, 1, 0)))
    c = random.uniform(10,50) * random_direction_close_to(matrix.col((0, 0, 1)))
    mycrystal = crystal(a, b, c)

    # make a dumb goniometer that rotates around X
    mygonio = goniometer_factory.known_axis((1, 0, 0))

    # generate some indices
    resolution = 1.0
    indices = full_sphere_indices(
        unit_cell = mycrystal.get_unit_cell(),
        resolution_limit = resolution,
        space_group = space_group(space_group_symbols(1).hall()))

    # generate list of phi values
    R_to_rossmann = align_reference_frame(
        mybeam.get_direction(), (0.0, 0.0, 1.0),
        mygonio.get_rotation_axis(), (0.0, 1.0, 0.0))

    ra = rotation_angles(resolution,
                    R_to_rossmann * mycrystal.get_U() * mycrystal.get_B(),
                    mybeam.get_wavelength(),
                    R_to_rossmann * matrix.col(mygonio.get_rotation_axis()))

    obs_indices, obs_angles = ra.observed_indices_and_angles_from_angle_range(
        phi_start_rad = 0.0, phi_end_rad = math.pi, indices = indices)

    # test rotation_angles wrapper
    ap = angle_predictor(mycrystal, mybeam, mygonio, resolution)
    obs_indices2, obs_angles2 = ap.observed_indices_and_angles_from_angle_range(
        phi_start_rad = 0.0, phi_end_rad = math.pi, indices = indices)

    # NB
    # * obs_indices is a scitbx_array_family_flex_ext.vec3_double and contains
    # floating point indices.
    # * obs_indices2 is a cctbx_array_family_flex_ext.miller_index containing
    # integer indices

    for h1, h2, p1, p2 in zip(obs_indices, obs_indices2, obs_angles, obs_angles2):
        assert(h1 == h2)
        assert(p1 == p2)

    # Now we're ready to test angle_predictor_py
    ap = angle_predictor_py(mycrystal, mybeam, mygonio, resolution)

    for hkl, phi in zip(obs_indices2, obs_angles2):
        ap_phi = ap.predict(hkl)
        if ap_phi is None: continue
        assert (approx_equal(ap_phi[0], phi, out = None)) or \
               (approx_equal(ap_phi[1], phi, out = None))

    # Test prediction of multiple indices from angle_predictor_py
    obs_indices3, obs_angles3 = ap.observed_indices_and_angles_from_angle_range(
        phi_start_rad = 0.0, phi_end_rad = math.pi, indices = indices)

    # Test method for multiple indices
    # NB I'm going to do something a bit horrible here. The values in obs_indices
    # are floats (it is an array of vec3), while obs_indices2 are ints (a Miller
    # index array). To compare the sets of indices I convert obs_indices back to
    # integer Miller indices
    temp = flex.miller_index()
    for hkl in obs_indices:
        hkl = map(lambda x: int(round(x)), hkl)
        temp.append(tuple(hkl))
    obs_indices = temp

    for i in xrange(20):
        print (obs_indices[i], obs_angles[i], obs_indices2[i], obs_angles2[i],
              obs_indices3[i], obs_angles3[i])
    print len(obs_indices), len(obs_indices2), len(obs_indices3)

    # FIXME how do I compare these sets?
    # There are some differences between the limits of reflections that are
    # predicted by my Python reflection prediction code in angle_predictor_py
    # and those by the rotation_angles wrapper, angle_predictor. That can be
    # seen by running this script a few times and looking at the output lines
    # in the above loop. Sometimes one method predicts an extra 'shell' of
    # reflections compared to the other.

#    obs_indices2, obs_angles2 = ap.observed_indices_and_angles_from_angle_range(
#        phi_start_rad = 0.0, phi_end_rad = math.pi, indices = obs_indices1)
#
#    print len(obs_indices2)
#    for i in range(10):
#        print obs_indices[i], obs_angles[i], obs_indices2[i], obs_angles2[i]
#    for ang, ang2 in zip(obs_angles, obs_angles2):
#        assert approx_equal(ang, ang2)
#

    # Now compare dials.algorithms.spot_prediction with the above methods. I
    # want to switch over to using DIALS spot prediction, but should verify that
    # it gives me the same results first.
    from dials.algorithms.spot_prediction import IndexGenerator, RotationAngles

    # Construct an IndexGenerator(unit_cell, space_group_type, d_min)
    index_generator = IndexGenerator(mycrystal.get_unit_cell(), space_group(space_group_symbols(1).hall()).type(), 1.0)
    dials_indices = index_generator.to_array()
    
    # Ensure generated indices are the same
    for h1, h2 in zip(indices, dials_indices):
        assert h1 == h2

    # Construct a RotationAngles(beam_vector, spindle_axis)
    rotation_angles = RotationAngles(mybeam.get_s0(),
                                     mygonio.get_rotation_axis())


    UB = mycrystal.get_U() * mycrystal.get_B()

    # try a few indices
    dials_angles = flex.double()
    dials_obs_indices = flex.miller_index()
    for h in indices:
        try:
            angles = rotation_angles(h, UB)
        except RuntimeError as e:
            continue
        dials_obs_indices.append(h)
        dials_angles.append(angles[0])
        dials_obs_indices.append(h)
        dials_angles.append(angles[1])

    for i in xrange(20):
        print (obs_indices2[i], obs_angles2[i], dials_obs_indices[i], dials_angles[i])
    print len(obs_indices2), len(dials_obs_indices)
    
    # Note dials_indices is about twice as long as obs_indices2, because it
    # hasn't been filtered by a rotation range of pi radians!

    # if we got this far,
    print "OK"
