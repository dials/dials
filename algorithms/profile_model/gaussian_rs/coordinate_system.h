/*
 * coordinate_system.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_PROFILE_MODEL_GAUSSIAN_RS_COORDINATE_SYSTEM_H
#define DIALS_ALGORITHMS_PROFILE_MODEL_GAUSSIAN_RS_COORDINATE_SYSTEM_H

#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/array_family/tiny_types.h>
#include <dxtbx/model/panel.h>
#include <dials/error.h>

namespace dials {
  namespace algorithms {
    namespace profile_model {
      namespace gaussian_rs {

  using dxtbx::model::plane_ray_intersection;
  using dxtbx::model::plane_world_coordinate;
  using scitbx::mat3;
  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::af::double4;

  /**
   * Helper function to calculate path length correction factor.
   * @param m2 The rotation axis
   * @param e1 The e1 axis of the reflection coordinate system
   * @returns Zeta the path length correction factor.
   */
  inline double zeta_factor(vec3<double> m2, vec3<double> e1) {
    return m2 * e1;
  }

  /**
   * Helper function to calculate path length correction factor.
   * @param m2 The rotation axis
   * @param s0 The incident beam vector
   * @param s1 The diffracted beam vector
   * @returns Zeta the path length correction factor.
   */
  inline double zeta_factor(vec3<double> m2, vec3<double> s0, vec3<double> s1) {
    vec3<double> e1 = s1.cross(s0);
    DIALS_ASSERT(e1.length() > 0);
    return zeta_factor(m2, e1.normalize());
  }

  /**
   * Class representing the local reflection coordinate system
   */
  class CoordinateSystem2d {
  public:
    /**
     * Initialise coordinate system. s0 should be the same length as s1.
     * These quantities are not checked because this class will be created for
     * each reflection and we want to maximize performance.
     * @param m2 The rotation axis
     * @param s0 The incident beam vector
     * @param s1 The diffracted beam vector
     * @param phi The rotation angle
     */
    CoordinateSystem2d(vec3<double> s0, vec3<double> s1)
        : s0_(s0),
          s1_(s1.normalize() * s0.length()),
          p_star_(s1 - s0),
          e1_(s1.cross(s0).normalize()),
          e2_(s1.cross(e1_).normalize()) {}

    /** @returns the incident beam vector */
    vec3<double> s0() const {
      return s0_;
    }

    /** @returns the diffracted beam vector */
    vec3<double> s1() const {
      return s1_;
    }

    /** @returns the rotated reciprocal space vector */
    vec3<double> p_star() const {
      return p_star_;
    }

    /** @returns the e1 axis vector */
    vec3<double> e1_axis() const {
      return e1_;
    }

    /** @returns the e2 axis vector */
    vec3<double> e2_axis() const {
      return e2_;
    }

    /**
     * Transform the beam vector to the reciprocal space coordinate system.
     * @param s_dash The beam vector
     * @returns The e1, e2 coordinates
     */
    vec2<double> from_beam_vector(const vec3<double> &s_dash) const {
      double s1_length = s1_.length();
      DIALS_ASSERT(s1_length > 0);
      vec3<double> scaled_e1 = e1_ / s1_length;
      vec3<double> scaled_e2 = e2_ / s1_length;
      return vec2<double>(scaled_e1 * (s_dash - s1_), scaled_e2 * (s_dash - s1_));
    }

    /**
     * Transform the reciprocal space coordinate to get the beam vector.
     * @param c12 The e1 and e2 coordinates.
     * @returns The beam vector
     */
    vec3<double> to_beam_vector(const vec2<double> &c12) const {
      double radius = s1_.length();
      DIALS_ASSERT(radius > 0);
      vec3<double> scaled_e1 = e1_ * radius;
      vec3<double> scaled_e2 = e2_ * radius;
      vec3<double> normalized_s1 = s1_ / radius;

      vec3<double> p = c12[0] * scaled_e1 + c12[1] * scaled_e2;
      double b = radius * radius - p.length_sq();
      DIALS_ASSERT(b >= 0);
      double d = -(normalized_s1 * p) + std::sqrt(b);
      return p + d * normalized_s1;
    }

  private:
    vec3<double> s0_;
    vec3<double> s1_;
    vec3<double> p_star_;
    vec3<double> e1_;
    vec3<double> e2_;
  };

  /**
   * Class representing the local reflection coordinate system
   */
  class CoordinateSystem {
  public:
    /**
     * Initialise coordinate system. s0 should be the same length as s1.
     * These quantities are not checked because this class will be created for
     * each reflection and we want to maximize performance.
     * @param m2 The rotation axis
     * @param s0 The incident beam vector
     * @param s1 The diffracted beam vector
     * @param phi The rotation angle
     */
    CoordinateSystem(vec3<double> m2, vec3<double> s0, vec3<double> s1, double phi)
        : m2_(m2.normalize()),
          s0_(s0),
          s1_(s1),
          phi_(phi),
          p_star_(s1 - s0),
          e1_(s1.cross(s0).normalize()),
          e2_(s1.cross(e1_).normalize()),
          e3_((s1 + s0).normalize()),
          zeta_(zeta_factor(m2_, e1_)) {}

    /** @returns The rotation axis */
    vec3<double> m2() const {
      return m2_;
    }

    /** @returns the incident beam vector */
    vec3<double> s0() const {
      return s0_;
    }

    /** @returns the diffracted beam vector */
    vec3<double> s1() const {
      return s1_;
    }

    /** @returns the rotation angle */
    double phi() const {
      return phi_;
    }

    /** @returns the rotated reciprocal space vector */
    vec3<double> p_star() const {
      return p_star_;
    }

    /** @returns the e1 axis vector */
    vec3<double> e1_axis() const {
      return e1_;
    }

    /** @returns the e2 axis vector */
    vec3<double> e2_axis() const {
      return e2_;
    }

    /** @returns the e3 axis vector */
    vec3<double> e3_axis() const {
      return e3_;
    }

    /** @returns the zeta factor (related to lorentz correction) */
    double zeta() const {
      return zeta_;
    }

    /** @returns the inverse lorentz correction factor */
    double lorentz_inv() const {
      return std::abs(m2_ * (s1_.cross(s0_))) / (s0_.length() * s1_.length());
    }

    /** @returns the lorentz correction factor */
    double lorentz() const {
      return 1.0 / lorentz_inv();
    }

    /** @returns the increase in the length of the shortest path */
    double path_length_increase() const {
      DIALS_ASSERT(zeta_ != 0.0);
      return std::abs(1.0 / zeta_);
    }

    /**
     * Return the limits of the acceptable coordinates in the reciprocal
     * space coordinate system. The values are returned as a pair in the
     * following form.
     *  (sqrt(max_e1_length**2 + max_e2_length**2), max_e3_length)
     *
     * The e3 limit is found by rotating the p* vector by pi/2 and projecting
     * onto the e3 axis. This equation can be simplified to the one below
     * which only requires the projection of the rotation axis onto the
     * e1, e3 and p* vectors giving: m2xe1 + (m2xe3)(m2xp*)|p*|
     *
     * @returns The limits
     */
    double4 limits() const {
      vec3<double> p_star_norm = p_star_.normalize();
      double m2e3 = m2_ * e3_;
      double m2ps = m2_ * p_star_norm;
      double m2e1 = m2_ * e1_;
      double m2e1_m2e1 = m2e1 * m2e1;
      double m2e3_m2ps = m2e3 * m2ps;
      double r = m2e3_m2ps * m2e3_m2ps + m2e1_m2e1;
      DIALS_ASSERT(r >= 0.0);
      return double4(-1.0, 1.0, m2e3_m2ps - std::sqrt(r), m2e3_m2ps + std::sqrt(r));
    }

    /**
     * Transform the beam vector to the reciprocal space coordinate system.
     * @param s_dash The beam vector
     * @returns The e1, e2 coordinates
     */
    vec2<double> from_beam_vector(const vec3<double> &s_dash) const {
      double s1_length = s1_.length();
      DIALS_ASSERT(s1_length > 0);
      vec3<double> scaled_e1 = e1_ / s1_length;
      vec3<double> scaled_e2 = e2_ / s1_length;
      return vec2<double>(scaled_e1 * (s_dash - s1_), scaled_e2 * (s_dash - s1_));
    }

    /**
     * Transform the rotation angle to the reciprocal space coordinate system
     * @param phi_dash The rotation angle
     * @returns The e3 coordinate.
     */
    double from_rotation_angle(double phi_dash) const {
      double p_star_length = p_star_.length();
      DIALS_ASSERT(p_star_length > 0);
      vec3<double> scaled_e3 = e3_ / p_star_length;
      return scaled_e3
             * (p_star_.unit_rotate_around_origin(m2_, phi_dash - phi_) - p_star_);
    }

    /**
     * Transform the rotation to the reciprocal space coordinate system using
     * a fast approximate method.
     * @param phi_dash The rotation angle
     * @retuns The e3 coordinate
     */
    double from_rotation_angle_fast(double phi_dash) const {
      return zeta_ * (phi_dash - phi_);
    }

    /**
     * Transform the beam vector and rotation angle to get the full
     * reciprocal space coordinate
     * @param s1_dash The beam vector
     * @param phi_dash The rotation angle
     * @returns The reciprocal space coordinate
     */
    vec3<double> from_beam_vector_and_rotation_angle(vec3<double> s1_dash,
                                                     double phi_dash) const {
      vec2<double> c12 = from_beam_vector(s1_dash);
      return vec3<double>(c12[0], c12[1], from_rotation_angle_fast(phi_dash));
    }

    /**
     * Transform the reciprocal space coordinate to get the beam vector.
     * @param c12 The e1 and e2 coordinates.
     * @returns The beam vector
     */
    vec3<double> to_beam_vector(const vec2<double> &c12) const {
      double radius = s1_.length();
      DIALS_ASSERT(radius > 0);
      vec3<double> scaled_e1 = e1_ * radius;
      vec3<double> scaled_e2 = e2_ * radius;
      vec3<double> normalized_s1 = s1_ / radius;

      vec3<double> p = c12[0] * scaled_e1 + c12[1] * scaled_e2;
      double b = radius * radius - p.length_sq();
      DIALS_ASSERT(b >= 0);
      double d = -(normalized_s1 * p) + std::sqrt(b);
      return p + d * normalized_s1;
    }

    /**
     * Apply the transform by solving the following equation for t
     *  c3 = (m2.e1)sin(dt) + (m2.e3)(m2.p*)(1 - cos(dt))
     * Giving:
     *  dt = 2 atan((sqrt((m2.e1)^2 + 2 c3(m2.e3)(m2.p*) - c3^2) + m2.e1) /
     *             c3 - 2 (m2.e3)(m2.p*))
     *
     * @param c3 The XDS e3 coordinate
     * @returns The rotation angle phi'
     */
    double to_rotation_angle(double c3) const {
      double m2e1 = m2_ * e1_;
      double m2e1_m2e1 = (m2e1 * m2e1);
      double m2e3_m2ps = (2.0 * (m2_ * e3_) * (m2_ * p_star_.normalize()));

      double l = m2e1_m2e1 + c3 * m2e3_m2ps - c3 * c3;
      DIALS_ASSERT(l >= 0.0);
      double y = std::sqrt(l) + m2e1;
      double x = c3 - m2e3_m2ps;
      DIALS_ASSERT(x != 0.0);
      return phi_ + 2.0 * atan(y / x);
    }

    /**
     * Transform the e3 coordinate to get the rotation angle using a
     * fast approximate method.
     * @param c3 The e3 coordinate
     * @returns The rotation angle.
     */
    double to_rotation_angle_fast(double c3) const {
      return c3 / zeta_ + phi_;
    }

    /**
     * Transform the reciprocal space coordinate back to real space to
     * get the beam vector and rotation angle.
     * @param The reciprocal space coordinate
     * @returns a pair of the beam vector and rotation angle
     */
    std::pair<vec3<double>, double> to_beam_vector_and_rotation_angle(
      vec3<double> c) const {
      return std::make_pair(to_beam_vector(vec2<double>(c[0], c[1])),
                            to_rotation_angle_fast(c[2]));
    }

  private:
    vec3<double> m2_;
    vec3<double> s0_;
    vec3<double> s1_;
    double phi_;
    vec3<double> p_star_;
    vec3<double> e1_;
    vec3<double> e2_;
    vec3<double> e3_;
    double zeta_;
  };

}}}}  // namespace dials::algorithms::profile_model::gaussian_rs

#endif  // DIALS_ALGORITHMS_PROFILE_MODEL_GAUSSIAN_RS_COORDINATE_SYSTEM_H
