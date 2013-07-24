/*
 * filter.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_INTEGRATION_FILTER_H
#define DIALS_ALGORITHMS_INTEGRATION_FILTER_H

#include <cmath>
#include <scitbx/vec3.h>
#include <dxtbx/model/goniometer.h>
#include <dxtbx/model/beam.h>
#include <dials/algorithms/reflection_basis/coordinate_system.h>
#include <dials/model/data/reflection.h>

namespace dials { namespace algorithms { namespace filter {

  using scitbx::vec3;
  using dxtbx::model::Goniometer;
  using dxtbx::model::Beam;
  using dials::algorithms::reflection_basis::CoordinateSystem;
  using dials::algorithms::reflection_basis::zeta_factor;
  using dials::model::Reflection;
  using dials::model::ReflectionList;

  /**
   * Calculate the zeta factor and check its absolute value is above the
   * minimum specified value.
   * @param m2 The rotation axis (normalized)
   * @param s0 The incident beam vector
   * @param s1 The diffracted beam vector
   * @param zeta_min The minimum allowed zeta value
   * @returns True/False, zeta is valid
   */
  inline
  bool is_zeta_valid(vec3<double> m2, vec3<double> s0, vec3<double> s1,
      double zeta_min) {
    return std::abs(zeta_factor(m2, s0, s1)) >= zeta_min;
  }

  /**
   * Calculate the zeta factor and check its absolute value is above the
   * minimum specified value.
   * @param cs The local reflection coordinate system
   * @param zeta_min The minimum allowed zeta value
   * @returns True/False, zeta is valid
   */
  inline
  bool is_zeta_valid(const CoordinateSystem &cs, double zeta_min) {
    return std::abs(cs.zeta()) >= zeta_min;
  }

  /**
   * Calculate the zeta factor and check its absolute value is above the
   * minimum specified value.
   * @param m2 The rotation axis (normalized)
   * @param s0 The incident beam vector
   * @param r The reflection
   * @param zeta_min The minimum allowed zeta value
   * @returns True/False, zeta is valid
   */
  inline
  bool is_zeta_valid(vec3<double> m2, vec3<double> s0, const Reflection &r,
      double zeta_min) {
    return is_zeta_valid(m2, s0, r.get_beam_vector(), zeta_min);
  }

  /**
   * Calculate the zeta factor and check its absolute value is above the
   * minimum specified value.
   * @param g The goniometer
   * @param b The beam
   * @param r The reflection
   * @param zeta_min The minimum allowed zeta value
   * @returns True/False, zeta is valid
   */
  inline
  bool is_zeta_valid(const Goniometer &g, const Beam &b, const Reflection &r,
      double zeta_min) {
    return is_zeta_valid(g.get_rotation_axis(), b.get_s0(), r, zeta_min);
  }

  /**
   * Check if the XDS small angle approximation holds for the local
   * reflection transform.
   *
   * Check that the following condition holds:
   *  (m2.e1)^2 + 2*c3*(m2.e3)*(m2.p*) - c3^2 >= 0
   *
   * @param m2 The rotation axis
   * @param s0 The incident beam vector
   * @param s1 The diffracted beam vector
   * @param delta_m The mosaicity * n_sigma
   * @returns True/False, the small angle approximation is valid
   */
  inline
  bool is_xds_small_angle_valid(vec3<double> m2, vec3<double> s0,
      vec3<double> s1, double delta_m) {
    vec3<double> ps = (s1 - s0).normalize();
    vec3<double> e1 = s1.cross(s0).normalize();
    vec3<double> e3 = (s1 + s0).normalize();
    double m2e1 = m2 * e1;
    double m2e3 = m2 * e3;
    double m2ps = m2 * ps;
    double c3 = -std::abs(delta_m);
    return (m2e1*m2e1 + 2.0*c3*m2e3*m2ps - c3*c3) >= 0.0;
  }

  /**
   * Check if the XDS small angle approximation holds for the local
   * reflection transform.
   *
   * Check that the following condition holds:
   *  (m2.e1)^2 + 2*c3*(m2.e3)*(m2.p*) - c3^2 >= 0
   *
   * @param cs The local reflection coordinate system
   * @param delta_m The mosaicity * n_sigma
   * @returns True/False, the small angle approximation is valid
   */
  inline
  bool is_xds_small_angle_valid(const CoordinateSystem &cs, double delta_m) {
    vec3<double> m2 = cs.m2();
    vec3<double> ps = cs.p_star().normalize();
    vec3<double> e1 = cs.e1_axis();
    vec3<double> e3 = cs.e3_axis();
    double m2e1 = m2 * e1;
    double m2e3 = m2 * e3;
    double m2ps = m2 * ps;
    double c3 = -std::abs(delta_m);
    return (m2e1*m2e1 + 2.0*c3*m2e3*m2ps - c3*c3) >= 0.0;
  }

  /**
   * Check if the XDS small angle approximation holds for the local
   * reflection transform.
   *
   * Check that the following condition holds:
   *  (m2.e1)^2 + 2*c3*(m2.e3)*(m2.p*) - c3^2 >= 0
   *
   * @param m2 The rotation axis
   * @param s0 The incident beam vector
   * @param r The reflection
   * @param delta_m The mosaicity * n_sigma
   * @returns True/False, the small angle approximation is valid
   */
  inline
  bool is_xds_small_angle_valid(vec3<double> m2, vec3<double> s0,
      const Reflection &r, double delta_m) {
    return is_xds_small_angle_valid(m2, s0, r.get_beam_vector(), delta_m);
  }

  /**
   * Check if the XDS small angle approximation holds for the local
   * reflection transform.
   *
   * Check that the following condition holds:
   *  (m2.e1)^2 + 2*c3*(m2.e3)*(m2.p*) - c3^2 >= 0
   *
   * @param g The goniometer
   * @param b The beam
   * @param r The reflection
   * @param delta_m The mosaicity * n_sigma
   * @returns True/False, the small angle approximation is valid
   */
  inline
  bool is_xds_small_angle_valid(const Goniometer &g, const Beam &b,
      const Reflection &r, double delta_m) {
    return is_xds_small_angle_valid(g.get_rotation_axis(), b.get_s0(),
      r, delta_m);
  }

  /**
   * Check that the angle can be mapped to the local reflection coordinate
   * system
   * @param m2 The rotation axis
   * @param s0 The incident beam vector
   * @param s1 The diffracted beam vector
   * @param delta_m The mosaicity * n_sigma
   * @returns True/False, the angle is valid
   */
  inline
  bool is_xds_angle_valid(vec3<double> m2, vec3<double> s0, vec3<double> s1,
      double delta_m) {
    vec3<double> ps = (s1 - s0).normalize();
    vec3<double> e1 = s1.cross(s0).normalize();
    vec3<double> e3 = (s1 + s0).normalize();
    double m2e1 = m2 * e1;
    double m2e3 = m2 * e3;
    double m2ps = m2 * ps;
    double m2e3_m2ps = m2e3 * m2ps;
    if (m2e1 == 0) {
      return false;
    }
    double rt = std::sqrt(m2e1*m2e1 + m2e3_m2ps*m2e3_m2ps);
    double tanphi0 = (m2e3_m2ps + rt) / m2e1;
    double tanphi1 = (m2e3_m2ps - rt) / m2e1;
    double phi0 = 2.0 * std::atan(tanphi0);
    double phi1 = 2.0 * std::atan(tanphi1);
    return phi0 <= delta_m && phi1 >= delta_m;
  }

  /**
   * Check that the angle can be mapped to the local reflection coordinate
   * system
   * @param cs The local reflection coordinate system
   * @param delta_m The mosaicity * n_sigma
   * @returns True/False, the angle is valid
   */
  inline
  bool is_xds_angle_valid(const CoordinateSystem &cs, double delta_m) {
    vec3<double> m2 = cs.m2();
    vec3<double> ps = cs.p_star().normalize();
    vec3<double> e1 = cs.e1_axis();
    vec3<double> e3 = cs.e3_axis();
    double m2e1 = m2 * e1;
    double m2e3 = m2 * e3;
    double m2ps = m2 * ps;
    double m2e3_m2ps = m2e3 * m2ps;
    if (m2e1 == 0) {
      return false;
    }
    double rt = std::sqrt(m2e1*m2e1 + m2e3_m2ps*m2e3_m2ps);
    double tanphi0 = (m2e3_m2ps + rt) / m2e1;
    double tanphi1 = (m2e3_m2ps - rt) / m2e1;
    double phi0 = 2.0 * std::atan(tanphi0);
    double phi1 = 2.0 * std::atan(tanphi1);
    return phi0 <= delta_m && phi1 >= delta_m;
  }

  /**
   * Check that the angle can be mapped to the local reflection coordinate
   * system
   * @param m2 The rotation axis
   * @param s0 The incident beam vector
   * @param r The reflection
   * @param delta_m The mosaicity * n_sigma
   * @returns True/False, the angle is valid
   */
  inline
  bool is_xds_angle_valid(vec3<double> m2, vec3<double> s0,
      const Reflection &r, double delta_m) {
    return is_xds_angle_valid(m2, s0, r.get_beam_vector(), delta_m);
  }

  /**
   * Check that the angle can be mapped to the local reflection coordinate
   * system
   * @param g The goniometer
   * @param s0 The beam
   * @param r The reflection
   * @param delta_m The mosaicity * n_sigma
   * @returns True/False, the angle is valid
   */
  inline
  bool is_xds_angle_valid(const Goniometer &g, const Beam &b,
      const Reflection &r, double delta_m) {
    return is_xds_angle_valid(g.get_rotation_axis(), b.get_s0(), r, delta_m);
  }

  /**
   * Filter the reflection list by the value of zeta. Set any reflections
   * below the value to invalid.
   * @param g The goniometer
   * @param b The beam
   * @param r The reflection list
   * @param min_zeta The minimum zeta value
   */
  void by_zeta(const Goniometer &g, const Beam &b, ReflectionList &r,
      double min_zeta) {
    for (std::size_t i = 0; i < r.size(); ++i) {
      if (!is_zeta_valid(g, b, r[i], min_zeta)) {
        r[i].set_valid(false);
      }
    }
  }

  /**
   * Filter the reflection list by the validity of the xds small angle approx.
   * Set any reflections for which its invalid to invalid.
   * @param g The goniometer
   * @param b The beam
   * @param r The reflection list
   * @param delta_m The mosaicity * n_sigma
   */
  void by_xds_small_angle(const Goniometer &g, const Beam &b,
      ReflectionList &r, double delta_m) {
    for (std::size_t i = 0; i < r.size(); ++i) {
      if (!is_xds_small_angle_valid(g, b, r[i], delta_m)) {
        r[i].set_valid(false);
      }
    }
  }

  /**
   * Filter the reflection list by the validity of the xds angle.
   * Set any reflections for which its invalid to invalid.
   * @param g The goniometer
   * @param b The beam
   * @param r The reflection list
   * @param delta_m The mosaicity * n_sigma
   */
  void by_xds_angle(const Goniometer &g, const Beam &b,
      ReflectionList &r, double delta_m) {
    for (std::size_t i = 0; i < r.size(); ++i) {
      if (!is_xds_angle_valid(g, b, r[i], delta_m)) {
        r[i].set_valid(false);
      }
    }
  }

}}} // namespace dials::algorithms::filter

#endif /* DIALS_ALGORITHMS_INTEGRATION_FILTER_H */
