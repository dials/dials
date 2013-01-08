
#ifndef DIALS_GEOMETRY_TRANSFORM_FROM_HKL_TO_BEAM_VECTOR_H
#define DIALS_GEOMETRY_TRANSFORM_FROM_HKL_TO_BEAM_VECTOR_H

#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include "../reciprocal_lattice_coordinate_system.h"

namespace dials { namespace geometry { namespace transform {

/** Class to represent geometry transform from miller indices to beam vector */
class from_hkl_to_beam_vector {

public:

    /** Default constructor */
    from_hkl_to_beam_vector() {}

    /** 
     * Initialise using the reciprocal lattice coordinate system 
     * @param rlcs The reciprocal lattice coordinate system class
     * @param s0 The incident beam vector
     * @param m2 The rotation axis
     */
    from_hkl_to_beam_vector(reciprocal_lattice_coordinate_system rlcs,
                            scitbx::vec3 <double> s0,
                            scitbx::vec3 <double> m2) 
        : _b1_star(rlcs.get_b1_star_axis()),
          _b2_star(rlcs.get_b2_star_axis()),
          _b3_star(rlcs.get_b3_star_axis()),
          _s0(s0),
          _m2(m2.normalize()) {}

public:

    /**
     * Apply the transform to a (h k l) point with a rotation
     * @param hkl The miller indices
     * @param phi The rotation angle
     * @returns The diffracted beam vector
     */
    scitbx::vec3 <double> apply(scitbx::vec3 <int> hkl, double phi) {
        return _s0 + (double(hkl[0]) * _b1_star + 
                      double(hkl[1]) * _b2_star + 
                      double(hkl[2]) * _b3_star)
                        .unit_rotate_around_origin(_m2, phi);
    }

private:

    scitbx::vec3 <double> _b1_star;
    scitbx::vec3 <double> _b2_star;
    scitbx::vec3 <double> _b3_star;
    scitbx::vec3 <double> _s0;
    scitbx::vec3 <double> _m2;
};

}}}

#endif // DIALS_GEOMETRY_TRANSFORM_FROM_HKL_TO_BEAM_VECTOR_H
