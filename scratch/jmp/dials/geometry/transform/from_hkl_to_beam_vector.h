
#ifndef DIALS_GEOMETRY_TRANSFORM_FROM_HKL_TO_BEAM_VECTOR_H
#define DIALS_GEOMETRY_TRANSFORM_FROM_HKL_TO_BEAM_VECTOR_H

#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <scitbx/array_family/flex_types.h>
#include <cctbx/miller.h>
#include "../../error.h"

namespace dials { namespace geometry { namespace transform {

typedef cctbx::miller::index <> miller_index;
typedef scitbx::af::flex <scitbx::vec3 <double> >::type flex_vec3_double;
typedef scitbx::af::flex <miller_index>::type flex_miller_index;

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
    FromHklToBeamVector(scitbx::mat3 <double> ub_matrix,
                        scitbx::vec3 <double> s0,
                        scitbx::vec3 <double> m2) 
        : b1_star_(ub_matrix[0], ub_matrix[3], ub_matrix[6]),
          b2_star_(ub_matrix[1], ub_matrix[4], ub_matrix[7]),
          b3_star_(ub_matrix[2], ub_matrix[5], ub_matrix[8]),
          s0_(s0),
          m2_(m2.normalize()) {}

public:

    /**
     * Apply the transform to a (h k l) point with a rotation
     * @param hkl The miller indices
     * @param phi The rotation angle
     * @returns The diffracted beam vector
     */
    scitbx::vec3 <double> apply(miller_index hkl, double phi) {
        return s0_ + (double(hkl[0]) * b1_star_ + 
                      double(hkl[1]) * b2_star_ + 
                      double(hkl[2]) * b3_star_)
                        .unit_rotate_around_origin(m2_, phi);
    }

    /**
     * Apply the transform to an array of hkl and phis
     * @param hkl An array of miller indices
     * @param phi An array of rotation angles
     * @returns An array of beam vectors
     * @throws An exception if the sizes of hkl and phi are not equal
     */
    flex_vec3_double apply(const flex_miller_index &hkl, 
                           const scitbx::af::flex_double &phi) {
        DIALS_ASSERT(hkl.size() == phi.size());
        flex_vec3_double result(hkl.size());
        for (int i = 0; i < hkl.size(); ++i) {
            result[i] = apply(hkl[i], phi[i]);
        }
        return result;
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
