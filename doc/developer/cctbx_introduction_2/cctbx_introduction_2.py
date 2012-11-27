#!/usr/bin/env cctbx.python
# cctbx_introduction_2.py
#
# An illustration of how to use cctbx code with orientation matrices, and
# how to perform a simple geometric strategy calculation, thus demonstrating
# symmetry... and also showing how B matrices are not to be trusted...
#
# This will also be useful:
#                                                   /!
#                      Y-axis                      / !
#                        ^                        /  !
#                        !                       /   !
#                        !                      /    !
#                        !   /                 /  Xd !
#                        !  /                 / * ^  !
#                        ! /                  ! 3 !  !
#                        !/      X-ray beam   !   !  !
#                        /-----------------------/--!---->X-axis
#                       /                     !  / *1!
#                    <-/-                     ! /    !
#                     /  \+ve phi             ! Yd  /
#                    /   /                    ! 2  /
#                   /                         ! * /
#                  Z-axis                  Ys ^ _/ 
#                Rotation                     ! /| Xs
#                 axis                        !/
#                                             O          

def mosflm_to_rossmann(a_matrix):
    '''Convert A* matrix in Mosflm convention to Rossmann convention.'''

    from scitbx import matrix
    
    S = matrix.sqr((0, 0, 1, 1, 0, 0, 0, 1, 0))

    return S * a_matrix

def parse_mosflm_matrix(matrix_file):
    '''Parse the mosflm matrix file to get: U, B, wavelength, unit cell,
    returning A* = U B in the Rossmann coordinate frame.'''

    from scitbx import matrix
    from cctbx import uctbx

    tokens = open(matrix_file, 'r').read(512).replace('-', ' -').split()

    assert(len(tokens) in (30, 32))

    lUB = matrix.sqr(map(float, tokens[0:9]))
    U = matrix.sqr(map(float, tokens[12:21]))
    
    uc = uctbx.unit_cell(map(float, tokens[21:27]))

    # derive the wavelength

    A = lUB.inverse().elems
    a = matrix.col(A[0:3])
    b = matrix.col(A[3:6])
    c = matrix.col(A[6:9])

    wavelength = (uc.parameters()[0] / a.length() +
                  uc.parameters()[1] / b.length() +
                  uc.parameters()[2] / c.length()) / 3

    # now verify that we do not reproduce the lUB matrix - this is the
    # PDB definition...

    B = matrix.sqr(uc.fractionalization_matrix())

    # ... and this is the Mosflm version (implied)

    mB = (1.0 / wavelength) * U.inverse() * lUB

    # not important - let's just return A* (i.e. a*, b*, c*)

    return mosflm_to_rossmann((1.0 / wavelength) * lUB), uc, wavelength

def generate_reflection_indices(uc, dmin):
    '''Generate the possible reflection indices from a unit cell object: N.B.
    that these are not symmetry reduced.'''

    maxh, maxk, maxl = uc.max_miller_indices(dmin)

    indices = []

    for h in range(-maxh, maxh + 1):
        for k in range(-maxk, maxk + 1):
            for l in range(-maxl, maxl + 1):

                if h == 0 and k == 0 and l == 0:
                    continue

                if uc.d((h, k, l)) < dmin:
                    continue

                indices.append((h, k, l))

    return indices    

def remove_absent_indices(indices, space_group_number):
    '''From the given list of indices, remove those reflections which should
    be systematic absences according to the given space group.'''

    from cctbx.sgtbx import space_group, space_group_symbols

    sg = space_group(space_group_symbols(space_group_number).hall())

    present = []

    for hkl in indices:
        if not sg.is_sys_absent(hkl):
            present.append(hkl)

    return present

def generate_intersection_angles(a_matrix, dmin, wavelength, indices):
    '''From an A matrix following the Mosflm convention and the list of
    indices, return a list of phi, (h, k, l) where (typically) there will be
    two records corresponding to each h, k, l.'''

    from rstbx.diffraction import rotation_angles
    from scitbx import matrix
    import math

    ra = rotation_angles(dmin, a_matrix, wavelength, matrix.col((0, 1, 0)))

    phi_hkl = []

    r2d = 180.0 / math.pi

    for i in indices:
        if ra(i):
            phis = ra.get_intersection_angles()
            phi_hkl.append((phis[0] * r2d % 360, i))
            phi_hkl.append((phis[1] * r2d % 360, i))

    return phi_hkl

def select_reflections(phi0, phi1, phi_hkl):
    '''Select reflections in range phi0 to phi1 inclusive.'''

    return [ph[1] for ph in phi_hkl if (ph[0] >= phi0 and ph[0] <= phi1)]

def reduce_reflections_to_asu(space_group_number, indices):
    '''Reduce reflection indices to asymmetric unit.'''

    from cctbx.sgtbx import space_group, space_group_symbols
    from cctbx.array_family import flex
    from cctbx.miller import map_to_asu

    sg = space_group(space_group_symbols(space_group_number).hall())

    miller = flex.miller_index(indices)

    map_to_asu(sg.type(), False, miller)

    return [hkl for hkl in miller]

def strategy(a_matrix, dmin, symmetry):
    '''Compute which 45 degree wedge gives the best completeness of data.'''

    a_star, uc, wavelength = parse_mosflm_matrix(a_matrix)
    indices = generate_reflection_indices(uc, dmin)
    present = remove_absent_indices(indices, symmetry)
    n_unique = len(set(reduce_reflections_to_asu(symmetry, indices)))
    observable = generate_intersection_angles(a_star, dmin, wavelength,
                                              present)
    dphi = 45.0
    for phi0 in 0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0:
        observed = select_reflections(phi0, phi0 + dphi, observable)
        n_unique_wedge = len(set(reduce_reflections_to_asu(symmetry,
                                                           observed)))
        print '%6.2f %6.4f' % (phi0, float(n_unique_wedge) / float(n_unique))

    return
 

if __name__ == '__main__':

    import sys

    strategy(sys.argv[1], float(sys.argv[2]), int(sys.argv[3]))

