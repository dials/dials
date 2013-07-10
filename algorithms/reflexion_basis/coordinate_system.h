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

#include <scitbx/vec3.h>
#include <dials/error.h>

namespace dials { namespace algorithms { namespace reflexion_basis {

  using scitbx::vec3;

  /**
   * Helper function to calculate path length correction factor.
   * @param m2 The rotation axis
   * @param e1 The e1 axis of the reflexion coordinate system
   * @returns Zeta the path length correction factor.
   */
  inline
  double zeta_factor(vec3<double> m2, vec2<double> e1) {
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
      : m2_(m2),
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

    vec3<double> scaled_e1_;
    vec3<double> scaled_e2_;
    vec3<double> s1_;
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
      return zeta_ * (phi_dash - phi_))
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
        e3_(cs.e3_axis()),
        p_star_(cs.p_star()),
        phi_(cs.phi()) {}

    /**
     * Map the rotation angle to the e3 coordinate of the local basis
     * @param phi The angle.
     * @returns The e3 angle.
     */
    double operator()(double phi_dash) const {
      return scaled_e3_ * (p_star_.unit_rotate_around_origin(
        m2_, phi_ - phi_dash) - p_star_);
    }

  private:

    vec3<double> m2_;
    vec3<double> e3_;
    vec3<double> p_star_;
    double phi_;
  };


  /**
   * Class to transform a beam vector and rotation angle to the reciprocal
   * space coordinate system.
   * @tparam FromBeamVectorType The class to transform the beam vector
   * @tparam FromRotationAngleType The class to transform the rotation angle
   */
  template <typename FromBeamVectorType,
            typename FromRotationAngleType>
  class FromBeamVectorAndRotationAngle {
  public:

    // Useful typedefs
    typedef FromBeamVectorType from_beam_vector_type;
    typedef FromRotationAngleType from_rotation_angle_type;

    /**
     * Initialise the transform with the coordinate system
     * @param cs The coordinate system
     */
    FromBeamVectorAndRotationAngle(const CoordinateSystem &cs)
      : from_beam_vector_(cs),
        from_rotation_angle_(cs) {}

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
    from_beam_vector_type from_beam_vector_;
    from_rotation_angle_type from_rotation_angle_;
  };


  /** Mapping using fast algorithm */
  typedef FromBeamVectorAndRotationAngle<
    FromBeamVector, FromRotationAngleFast>
      FromBeamVectorAndRotationAngleFast;

  /** Mapping using accurate algorithm */
  typedef FromBeamVectorAndRotationAngle<
    FromBeamVector, FromRotationAngleAccurate>
      FromBeamVectorAndRotationAngleAccurate;


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
       phi_(cs.phi()) {} {}

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
      : radius_(cs.p_star().length()),
        phi_(cs.phi()) {}

    /**
     * Apply the transform by projecting the point onto the arc of rotation
     * of the reciprocal space vector p. Then calculate the angle
     * @param e3 The XDS e3 coordinate
     * @returns The rotation angle phi'
     */
    double operator()(double e3) const {
      double b = radius_ * radius_ - e3 * e3;
      DIALS_ASSERT(b >= 0.0);
      double y = radius_ - std::sqrt(b);
      return phi_ - std::atan2(y, e3);
    }

  private:
    double radius_;
    double phi_;
  };


  /**
   * Class to transform a the reciprocal space coordinate system point to
   * a beam vector and rotation angle
   * @tparam ToBeamVectorType The class to transform the beam vector
   * @tparam ToRotationAngleType The class to transform the rotation angle
   */
  template <typename ToBeamVectorType,
            typename ToRotationAngleType>
  class ToBeamVectorAndRotationAngle {
  public:

    // Useful typedefs
    typedef ToBeamVectorType to_beam_vector_type;
    typedef ToRotationAngleType to_rotation_angle_type;

    /**
     * Initialise the transform with the coordinate system
     * @param cs The coordinate system
     */
    ToBeamVectorAndRotationAngle(const CoordinateSystem &cs)
      : to_beam_vector_(cs),
        to_rotation_angle_(cs) {}

    /**
     * Map the point to rotation angle and beam vector
     * @param e The point
     * @returns The beam vector and rotation angle
     */
    std::pair<vec3<double>, double> operator()(vec3<double> e) const {
      return std::make_pair(
        to_beam_vector(vec2<double>(e[0], e[1])),
        to_rotation_angle(e[2]));
    }

  private:
    to_beam_vector_type to_beam_vector_;
    to_rotation_angle_type to_rotation_angle_;
  };


  /** Mapping using fast algorithm */
  typedef ToBeamVectorAndRotationAngle<
    ToBeamVector, ToRotationAngleFast>
      ToBeamVectorAndRotationAngleFast;

  /** Mapping using accurate algorithm */
  typedef ToBeamVectorAndRotationAngle<
    ToBeamVector, ToRotationAngleAccurate>
      ToBeamVectorAndRotationAngleAccurate;

}} // namespace = dials::algorithms::reflexion_basis

#endif // DIALS_ALGORITHMS_REFLEXION_BASIS_COORDINATE_SYSTEM_H
