
#ifndef DIALS_GEOMETRY_TRANSFORM_FROM_HKL_TO_BEAM_VECTOR_H
#define DIALS_GEOMETRY_TRANSFORM_FROM_HKL_TO_BEAM_VECTOR_H

#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include "../reciprocal_lattice_coordinate_system.h"

namespace dials { namespace geometry { namespace transform {

/** Class to represent geometry transform from miller indices to beam vector */
class FromHklToBeamVector {

public:

    /** Default constructor */
    FromHklToBeamVector() {}

    /** 
     * Initialise using the reciprocal lattice coordinate system 
     * @param rlcs The reciprocal lattice coordinate system class
     * @param s0 The incident beam vector
     * @param m2 The rotation axis
     */
    FromHklToBeamVector(ReciprocalLatticeCoordinateSystem rlcs,
                        scitbx::vec3 <double> s0,
                        scitbx::vec3 <double> m2) 
        : b1_star_(rlcs.get_b1_star_axis()),
          b2_star_(rlcs.get_b2_star_axis()),
          b3_star_(rlcs.get_b3_star_axis()),
          s0_(s0),
          m2_(m2.normalize()) {}

public:

    /**
     * Apply the transform to a (h k l) point with a rotation
     * @param hkl The miller indices
     * @param phi The rotation angle
     * @returns The diffracted beam vector
     */
    scitbx::vec3 <double> apply(scitbx::vec3 <int> hkl, double phi) {
        return s0_ + (double(hkl[0]) * b1_star_ + 
                      double(hkl[1]) * b2_star_ + 
                      double(hkl[2]) * b3_star_)
                        .unit_rotate_around_origin(m2_, phi);
    }

private:

    scitbx::vec3 <double> b1_star_;
    scitbx::vec3 <double> b2_star_;
    scitbx::vec3 <double> b3_star_;
    scitbx::vec3 <double> s0_;
    scitbx::vec3 <double> m2_;
};

}}} // namespace = dials::geometry::transform

#endif // DIALS_GEOMETRY_TRANSFORM_FROM_HKL_TO_BEAM_VECTOR_H
