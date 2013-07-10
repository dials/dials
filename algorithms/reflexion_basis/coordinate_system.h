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
#ifndef DIALS_ALGORITHMS_REFLEXION_BASIS_COORDINATE_SYSTEM_H
#define DIALS_ALGORITHMS_REFLEXION_BASIS_COORDINATE_SYSTEM_H

#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/array_family/tiny_types.h>
#include <dxtbx/model/panel.h>
#include <dials/error.h>

namespace dials { namespace algorithms { namespace reflexion_basis {

  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::mat3;
  using scitbx::af::double4;
  using dxtbx::model::plane_ray_intersection;
  using dxtbx::model::plane_world_coordinate;

  /**
   * Helper function to calculate path length correction factor.
   * @param m2 The rotation axis
   * @param e1 The e1 axis of the reflexion coordinate system
   * @returns Zeta the path length correction factor.
   */
  inline
  double zeta_factor(vec3<double> m2, vec3<double> e1) {
    return m2 * e1;
  }

  /**
   * Helper function to calculate path length correction factor.
   * @param m2 The rotation axis
   * @param s0 The incident beam vector
   * @param s1 The diffracted beam vector
   * @returns Zeta the path length correction factor.
   */
  inline
  double zeta_factor(vec3<double> m2, vec3<double> s0, vec3<double> s1) {
    return zeta_factor(m2, s1.cross(s0).normalize());
  }

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
    CoordinateSystem(vec3<double> m2, vec3<double> s0,
                     vec3<double> s1, double phi)
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

    /** @returns the lorentz correction factor */
    double zeta() const {
      return zeta_;
    }

    /** @returns the increase in the length of the shortest path */
    double path_length_increase() const {
      DIALS_ASSERT(zeta_ != 0.0);
      return abs(1.0 / zeta_);
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
      double m2ps = m2_ * p_star_.normalize();
      double m2e1 = m2_ * e1_;
      double m2e1_m2e1 = m2e1 * m2e1;
      double m2e3_m2ps = m2e3 * m2ps;
      double r = m2e3_m2ps * m2e3_m2ps + m2e1_m2e1;
      DIALS_ASSERT(r >= 0.0);
      return double4(
        -1.0,
         1.0,
        m2e3_m2ps - sqrt(r),
        m2e3_m2ps + sqrt(r));
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

  /**
   * Class to represent a geometry transform from beam vector to the local
   * reciprocal space coordinate system (e1 and e2).
   */
  class FromBeamVector {

  public:

    /**
     * Initialise the transform using the XDS coordinate system.
     * @param cs The XDS coordinate system
     * @param s1 The diffracted beam vector
     */
    FromBeamVector(const CoordinateSystem &cs)
      : s1_(cs.s1()),
        scaled_e1_(cs.e1_axis() / s1_.length()),
        scaled_e2_(cs.e2_axis() / s1_.length()) {}

    /**
     * Apply the transform to a beam vector
     * @param s_dash The diffracted beam vector to transform
     * @returns The point in XDS coordinates
     */
    vec2 <double> operator()(vec3<double> s_dash) const {
      return vec2 <double> (
        scaled_e1_ * (s_dash - s1_),
        scaled_e2_ * (s_dash - s1_));
    }

  private:

    vec3<double> s1_;
    vec3<double> scaled_e1_;
    vec3<double> scaled_e2_;
  };


  /**
   * Class to transform from rotation angle to the e3 coordinate of the
   * local reflexion coordinate system.
   */
  class FromRotationAngleFast {
  public:

    /**
     * Initialise the transform
     * @param cs The coordinate system.
     * @param m2 The rotation axis.
     * @param phi The rotation angle.
     */
    FromRotationAngleFast(const CoordinateSystem &cs)
      : phi_(cs.phi()),
        zeta_(cs.zeta()) {}

    /**
     * Map the rotation angle to the e3 coordinate of the local basis
     * @param phi The angle.
     * @returns The e3 angle.
     */
    double operator()(double phi_dash) const {
      return zeta_ * (phi_dash - phi_);
    }

  private:
    double phi_;
    double zeta_;
  };


  /**
   * Class to transform the rotation angle to the e3 coordinate using
   * a more accurate method.
   */
  class FromRotationAngleAccurate {
  public:

    /**
     * Initialise the class using the coordinate system
     * @param cs The reflextion coordinate system
     */
    FromRotationAngleAccurate(const CoordinateSystem &cs)
      : m2_(cs.m2().normalize()),
        scaled_e3_(cs.e3_axis() / cs.p_star().length()),
        p_star_(cs.p_star()),
        phi_(cs.phi()) {}

    /**
     * Map the rotation angle to the e3 coordinate of the local basis
     * @param phi The angle.
     * @returns The e3 angle.
     */
    double operator()(double phi_dash) const {
      return scaled_e3_ * (p_star_.unit_rotate_around_origin(
        m2_, phi_dash - phi_) - p_star_);
    }

  private:

    vec3<double> m2_;
    vec3<double> scaled_e3_;
    vec3<double> p_star_;
    double phi_;
  };


  /**
   * Class to transform a beam vector and rotation angle to the reciprocal
   * space coordinate system.
   */
  class FromBeamVectorAndRotationAngle {
  public:

    /**
     * Initialise the transform with the coordinate system
     * @param cs The coordinate system
     */
    FromBeamVectorAndRotationAngle(const CoordinateSystem &cs)
      : from_beam_vector_(cs),
        from_rotation_angle_(cs) {}

    /**
     * Map the rotation angle to the e3 coordinate of the local basis
     * @param phi The angle.
     * @returns The e3 angle.
     */
    double operator()(double phi_dash) const {
      return from_rotation_angle_(phi_dash);
    }

    /**
     * Apply the transform to a beam vector
     * @param s_dash The diffracted beam vector to transform
     * @returns The point in XDS coordinates
     */
    vec2 <double> operator()(vec3<double> s_dash) const {
      return from_beam_vector_(s_dash);
    }

    /**
     * Map the rotation angle and beam vector to the coordinate system
     * @param s_dash The beam vector
     * @param phi_dash The rotation angle
     */
    vec3<double> operator()(vec3<double> s_dash, double phi_dash) const {
      vec2<double> e12 = from_beam_vector_(s_dash);
      return vec3<double>(e12[0], e12[1], from_rotation_angle_(phi_dash));
    }

  private:
    FromBeamVector from_beam_vector_;
    FromRotationAngleFast from_rotation_angle_;
  };


  /**
   * A class to perform a transform from XDS coordinate frame to beam vector
   */
  class ToBeamVector {
  public:

    /**
     * Initialise the transform
     * @param cs The XDS coordinate system
     */
    ToBeamVector(const CoordinateSystem &cs)
      : radius_(cs.s1().length()),
        scaled_e1_(cs.e1_axis() * radius_),
        scaled_e2_(cs.e2_axis() * radius_),
        normalized_s1_(cs.s1() / radius_) {}

    /**
     * Apply the transform to the xds point. The transform is done by
     * calculating the coordinate of point in the sphere tangent plane at the
     * location of the exit point of the beam vector, s1. The beam vector, s',
     * is then calculated by finding the intersection of line defined by the
     * plane normal (i.e. s1 vector) with origin at the calculated tangent plane
     * point point and the ewald sphere.
     *
     * @param e12 The XDS coordinate
     * @returns The beam vector
     */
    vec3<double> operator()(vec2<double> e12) const {
      vec3 <double> p = e12[0] * scaled_e1_ + e12[1] * scaled_e2_;
      double b = radius_ * radius_ - p.length_sq();
      DIALS_ASSERT(b >= 0);
      double d = -(normalized_s1_ * p) + std::sqrt(b);
      return p + d * normalized_s1_;
    }

  private:
    double radius_;
    vec3 <double> scaled_e1_;
    vec3 <double> scaled_e2_;
    vec3 <double> normalized_s1_;
  };


  /**
   * A class to transform from XDS e3 coord to the rotation angle, phi
   */
  class ToRotationAngleFast {
  public:

    /**
     * Initialise the class
     * @param cs The coordinate system
     */
    ToRotationAngleFast(const CoordinateSystem &cs)
     : zeta_(cs.zeta()),
       phi_(cs.phi()) {}

    /**
     * Apply the transform
     * @param e3 The XDS e3 coordinate
     * @returns The rotation angle phi'
     */
    double operator()(double e3) const {
      return e3 / zeta_ + phi_;
    }

  private:
    double zeta_;
    double phi_;
  };


  /**
   * Class to do a more accurate mapping from e3 to phi.
   */
  class ToRotationAngleAccurate {
  public:

    /**
     * Initialise the class
     * @param cs The coordinate system
     */
    ToRotationAngleAccurate(const CoordinateSystem &cs)
      : m2e1_(cs.m2() * cs.e1_axis()),
        m2e1_m2e1_(m2e1_ * m2e1_),
        m2e3_m2ps_(2.0 * (cs.m2() * cs.e3_axis()) *
          (cs.m2() * cs.p_star().normalize())) {}

    /**
     * Apply the transform by solving the following equation for t
     *  e3 = (m2.e1)sin(t) + (m2.e3)(m2.p*)(1 - cos(t))
     * Giving:
     *  t = 2 atan((sqrt((m2.e1)^2 + 2 c3(m2.e3)(m2.p*) - c3^2) + m2.e1) /
     *             c3 - 2 (m2.e3)(m2.p*))
     *
     * @param e3 The XDS e3 coordinate
     * @returns The rotation angle phi'
     */
    double operator()(double c3) const {
      double l = m2e1_m2e1_ + c3 * m2e3_m2ps_ - c3*c3;
      DIALS_ASSERT(l >= 0.0);
      double y = sqrt(l) + m2e1_;
      double x = c3 - m2e3_m2ps_;
      DIALS_ASSERT(x != 0.0);
      return 2.0*atan(y / x);
    }

  private:
    double m2e1_;
    double m2e1_m2e1_;
    double m2e3_m2ps_;
  };


  /**
   * Class to transform a the reciprocal space coordinate system point to
   * a beam vector and rotation angle
   */
  class ToBeamVectorAndRotationAngle {
  public:

    /**
     * Initialise the transform with the coordinate system
     * @param cs The coordinate system
     */
    ToBeamVectorAndRotationAngle(const CoordinateSystem &cs)
      : to_beam_vector_(cs),
        to_rotation_angle_(cs) {}

    /**
     * Apply the transform by projecting the point onto the arc of rotation
     * of the reciprocal space vector p. Then calculate the angle
     * @param e3 The XDS e3 coordinate
     * @returns The rotation angle phi'
     */
    double operator()(double e3) const {
      return to_rotation_angle_(e3);
    }

    /**
     * Apply the transform to the xds point. The transform is done by
     * calculating the coordinate of point in the sphere tangent plane at the
     * location of the exit point of the beam vector, s1. The beam vector, s',
     * is then calculated by finding the intersection of line defined by the
     * plane normal (i.e. s1 vector) with origin at the calculated tangent plane
     * point point and the ewald sphere.
     *
     * @param e12 The XDS coordinate
     * @returns The beam vector
     */
    vec3<double> operator()(vec2<double> e12) const {
      return to_beam_vector_(e12);
    }

    /**
     * Map the point to rotation angle and beam vector
     * @param e The point
     * @returns The beam vector and rotation angle
     */
    std::pair<vec3<double>, double> operator()(vec3<double> e) const {
      return std::make_pair(
        to_beam_vector_(vec2<double>(e[0], e[1])),
        to_rotation_angle_(e[2]));
    }

  private:
    ToBeamVector to_beam_vector_;
    ToRotationAngleFast to_rotation_angle_;
  };


  /**
   * Class to map detector millimeter coordinates to the reciprocal space
   * coordinate system.
   */
  class FromDetector {
  public:

    /**
     * Initialise the transform with the coordinate system and the
     * detector d matrix
     * @param cs The coordinate system
     * @param d The detector d matrix
     */
    FromDetector(const CoordinateSystem &cs, mat3<double> d)
      : from_s1_phi_(cs),
        d_(d) {}

    /**
     * Map the detector millimeter coordinate to the e1/e2 coordinate.
     * @param xy The detector millimeter coordinate
     * @returns The e1/e2 coordinate in reciprocal space
     */
    vec2<double> operator()(vec2<double> xy) const {
      return from_s1_phi_(plane_world_coordinate(d_, xy));
    }

    /**
     * Map the detector coordinate and rotation angle to the coordinate system.
     * @param xy The detector millimeter coordinate
     * @param phi The rotation angle
     * @returns The local reciprocal space coordinate
     */
    vec3<double> operator()(vec2<double> xy, double phi) const {
      return from_s1_phi_(plane_world_coordinate(d_, xy), phi);
    }

  private:
    FromBeamVectorAndRotationAngle from_s1_phi_;
    mat3<double> d_;
  };


  /**
   * Class to map the reciprocal space coordinate system to detector
   * millimeter coordinates.
   */
  class ToDetector {
  public:

    /**
     * Initialise the transform with the coordinate system and the
     * detector D matrix
     * @param cs The coordinate system
     * @param D The detector D matrix
     */
    ToDetector(const CoordinateSystem &cs, mat3<double> D)
      : to_s1_phi_(cs),
        D_(D) {}

    /**
     * Map the local reciprocal space e1/e2 coordinate to detector mm.
     * @param e The e1/e2 coordinate
     * @returns The detector mm coordinate.
     */
    vec2<double> operator()(vec2<double> e) const {
      return plane_ray_intersection(D_, to_s1_phi_(e));
    }

    /**
     * Map the local reciprocal space coordinate to the detector millimeter
     * and rotation angle.
     * @param e The local reciprocal space coordinate
     * @returns The mm/rad detector, rotation angle.
     */
    std::pair<vec2<double>, double> operator()(vec3<double> e) const {
      std::pair<vec3<double>, double> s1_phi = to_s1_phi_(e);
      return std::make_pair(
        plane_ray_intersection(D_, s1_phi.first),
        s1_phi.second);
    }

  private:
    ToBeamVectorAndRotationAngle to_s1_phi_;
    mat3<double> D_;
  };

}}} // namespace = dials::algorithms::reflexion_basis

#endif // DIALS_ALGORITHMS_REFLEXION_BASIS_COORDINATE_SYSTEM_H
