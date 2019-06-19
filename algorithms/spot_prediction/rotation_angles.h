/*
 * rotation_angles.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_SPOT_PREDICTION_ROTATION_ANGLES_H
#define DIALS_ALGORITHMS_SPOT_PREDICTION_ROTATION_ANGLES_H

#include <scitbx/constants.h>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <cctbx/miller.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using scitbx::mat3;
  using scitbx::vec2;
  using scitbx::vec3;
  using std::atan2;
  using std::sqrt;

  /** Calculate the square of the value */
  template <typename T>
  T sqr(T &a) {
    return a * a;
  }

  /**
   * Calculate the rotation angles at which the point is in the diffracting
   * condition.
   */
  class RotationAngles {
  public:
    /**
     * Initialise the rotation angle calculator
     * @param s0 The incident beam vector
     * @param m2 The rotation axis
     */
    RotationAngles(vec3<double> s0, vec3<double> m2)
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
    vec2<double> operator()(vec3<double> pstar0) const {
      // Calculate sq length of pstar0 and ensure p*^2 <= 4s0^2
      double pstar0_len_sq = pstar0.length_sq();
      DIALS_ASSERT(pstar0_len_sq <= 4 * s0_.length_sq());

      // Calculate dot product of p*0 with m1 and m3
      double pstar0_d_m1 = pstar0 * m1_;
      double pstar0_d_m2 = pstar0 * m2_;
      double pstar0_d_m3 = pstar0 * m3_;

      // Calculate dot product of p* with m3
      double pstar_d_m3 = (-(0.5 * pstar0_len_sq) - (pstar0_d_m2 * s0_d_m2)) / s0_d_m3;

      // Calculate sq distance of p*0 from rotation axis and ensure that
      // rho^2 >= (p*.m3)^2
      double rho_sq = (pstar0_len_sq - sqr(pstar0_d_m2));
      DIALS_ASSERT(rho_sq >= sqr(pstar_d_m3));

      // Calculate dot product of p* with m1
      double pstar_d_m1 = sqrt(rho_sq - sqr(pstar_d_m3));

      // Calculate sin/cos of the two angles with +- p*.m1
      double sinphi1, sinphi2, cosphi1, cosphi2;
      cosphi1 = (+(pstar_d_m1 * pstar0_d_m1) + (pstar_d_m3 * pstar0_d_m3));
      cosphi2 = (-(pstar_d_m1 * pstar0_d_m1) + (pstar_d_m3 * pstar0_d_m3));
      sinphi1 = (+(pstar_d_m1 * pstar0_d_m3) - (pstar_d_m3 * pstar0_d_m1));
      sinphi2 = (-(pstar_d_m1 * pstar0_d_m3) - (pstar_d_m3 * pstar0_d_m1));

      // Return the two angles
      return vec2<double>(atan2(sinphi1, cosphi1), atan2(sinphi2, cosphi2));
    }

    /**
     * Helper function to calculate angles directly from miller index and
     * UB matrix.
     * @param miller_index The miller indices
     * @param ub_matrix The ub matrix
     * @returns The two rotation angles.
     */
    vec2<double> operator()(cctbx::miller::index<> miller_index,
                            mat3<double> ub_matrix) const {
      return operator()(ub_matrix *miller_index);
    }

  private:
    /** Calculate the goniometer m1 axis (m1 = m2xs0 / |m2xs0|) */
    vec3<double> calculate_goniometer_m1_axis() const {
      return m2_.cross(s0_).normalize();
    }

    /** Calculate the goniometer m3 axis (m3 = m1xs2) */
    vec3<double> calculate_goniometer_m3_axis() const {
      return m1_.cross(m2_).normalize();
    }

    vec3<double> s0_;
    vec3<double> m2_;
    vec3<double> m1_;
    vec3<double> m3_;
    double s0_d_m2;
    double s0_d_m3;
  };

  /**
   * Helper function
   */
  inline vec2<double> rotation_angles(vec3<double> s0,
                                      vec3<double> m2,
                                      vec3<double> r0) {
    RotationAngles rot(s0, m2);
    return rot(r0);
  }

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_SPOT_PREDICTION_ROTATION_ANGLES_H
