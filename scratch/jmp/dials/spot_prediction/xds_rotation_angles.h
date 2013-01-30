
#ifndef DIALS_SPOT_PREDICTION_XDS_ROTATION_ANGLES_H
#define DIALS_SPOT_PREDICTION_XDS_ROTATION_ANGLES_H

#include <scitbx/constants.h>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <scitbx/array_family/flex_types.h>
#include <cctbx/miller.h>
#include "../error.h"

namespace dials { namespace spot_prediction {

/** Calculate the square of the value */
template <typename T>
T sqr(T &a) {
    return a * a;
}

/**
 * Calculate the rotation angles at which the point is in the diffracting 
 * condition.
 */
class XdsRotationAngles {

public:

    /**
     * Initialise the rotation angle calculator
     * @param s0 The incident beam vector
     * @param m2 The rotation axis
     */
    XdsRotationAngles(scitbx::vec3 <double> s0, scitbx::vec3 <double> m2)
        : s0_(s0),
          m2_(m2.normalize()),
          m1_(calculate_goniometer_m1_axis()),
          m3_(calculate_goniometer_m3_axis()),
          s0_d_m2(s0_ * m2_),
          s0_d_m3(s0_ * m3_) {}
                              
    /**
     * Calculate the rotation angles using the XDS method
     * @param pstar0 The unrotated reciprocal space vector
     * @returns The two rotation angles that satisfy the laue equations
     * @throws error if no angles exist.
     */
    scitbx::vec2 <double> calculate(scitbx::vec3 <double> pstar0) const {

        // Calculate sq length of pstar0 and ensure p*^2 <= 4s0^2
        double pstar0_len_sq = pstar0.length_sq();
        DIALS_ASSERT(pstar0_len_sq <= 4 * s0_.length_sq());
        
        // Calculate dot product of p*0 with m1 and m3
        double pstar0_d_m1 = pstar0 * m1_;
        double pstar0_d_m2 = pstar0 * m2_;
        double pstar0_d_m3 = pstar0 * m3_;

        // Calculate dot product of p* with m3
        double pstar_d_m3 = (-(0.5 * pstar0_len_sq) -
                              (pstar0_d_m2 * s0_d_m2)) / s0_d_m3;

        // Calculate sq distance of p*0 from rotation axis and ensure that
        // rho^2 >= (p*.m3)^2
        double rho_sq = (pstar0_len_sq - sqr(pstar0_d_m2));
        DIALS_ASSERT(rho_sq >= sqr(pstar_d_m3));

        // Calculate dot product of p* with m1
        double pstar_d_m1 = std::sqrt(rho_sq - sqr(pstar_d_m3));

        // Calculate sin/cos of the two angles with +- p*.m1
        double sinphi1, sinphi2, cosphi1, cosphi2;
        cosphi1 = (+(pstar_d_m1 * pstar0_d_m1) + (pstar_d_m3 * pstar0_d_m3));
        cosphi2 = (-(pstar_d_m1 * pstar0_d_m1) + (pstar_d_m3 * pstar0_d_m3));
        sinphi1 = (+(pstar_d_m1 * pstar0_d_m3) - (pstar_d_m3 * pstar0_d_m1));
        sinphi2 = (-(pstar_d_m1 * pstar0_d_m3) - (pstar_d_m3 * pstar0_d_m1));

        // Return the two angles
        return scitbx::vec2 <double> (std::atan2(sinphi1, cosphi1), 
                                      std::atan2(sinphi2, cosphi2));
    }
    
    /**
     * Helper function to calculate angles directly from miller index and 
     * UB matrix.
     * @param miller_index The miller indices
     * @param ub_matrix The ub matrix
     * @returns The two rotation angles.
     */
    scitbx::vec2 <double> calculate(cctbx::miller::index <> miller_index, 
                                    scitbx::mat3 <double> ub_matrix) const {
        return calculate(ub_matrix * miller_index);
    }    

private:

    /** Calculate the goniometer m1 axis (m1 = m2xs0 / |m2xs0|) */
    scitbx::vec3 <double> calculate_goniometer_m1_axis() const {
        return m2_.cross(s0_).normalize();
    }

    /** Calculate the goniometer m3 axis (m3 = m1xs2) */
    scitbx::vec3 <double> calculate_goniometer_m3_axis() const {
        return m1_.cross(m2_).normalize();
    }

    scitbx::vec3 <double> s0_;    
    scitbx::vec3 <double> m2_;
    scitbx::vec3 <double> m1_;
    scitbx::vec3 <double> m3_;
    double s0_d_m2;
    double s0_d_m3;
};

/** Convert the angle mod 2PI */
inline 
double mod_2pi(double angle) {
    return angle - scitbx::constants::two_pi * 
        std::floor(angle / scitbx::constants::two_pi); 
}

/** Convert a pair of angles to mod 2PI */
inline 
scitbx::vec2 <double> mod_2pi(scitbx::vec2 <double> angles) {
	return scitbx::vec2 <double> (mod_2pi(angles[0]), mod_2pi(angles[1]));
}

/**
    * Check if the angle is within the filter range
    * @param angle The angle to check
    * @param range The angular range
    * @returns True/False the angle is in the filter range
    */
inline 
bool is_angle_in_range(double angle, scitbx::vec2 <double> range) {
    return mod_2pi(angle - range[1]) >= mod_2pi(angle - range[0])
        || mod_2pi(angle - range[1]) == 0
        || mod_2pi(range[0] - range[1]) == 0;
}

/**
 * Check if the array of angles are within the angular range.
 * @param angle The array of angles
 * @param range The angular range
 * @returns A boolean array that is true if the angle is in the range.
 */
inline
scitbx::af::flex_bool is_angle_in_range(const scitbx::af::flex_double &angle,
                                        scitbx::vec2 <double> range) {
    scitbx::af::flex_bool result(angle.size());
    for (int i = 0; i < angle.size(); ++i) {
        result[i] = is_angle_in_range(angle[i], range);
    }
    return result;
}

}} // namespace dials::spot_prediction

#endif // DIALS_SPOT_PREDICTION_XDS_ROTATION_ANGLES_H
