/*
 * xds_coordinate_system.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_INTEGRATION_XDS_COORDINATE_SYSTEM_H
#define DIALS_ALGORITHMS_INTEGRATION_XDS_COORDINATE_SYSTEM_H

#include <scitbx/vec3.h>

namespace dials { namespace algorithms {

  using scitbx::vec3;

  /**
   * Helper function to calculate path length correction factor.
   * @param s0 The incident beam vector
   * @param s1 The diffracted beam vector
   * @param m2 The rotation axis
   * @returns Zeta the path length correction factor.
   */
  inline
  double zeta_factor(vec3 <double> s0, vec3 <double> s1, vec3 <double> m2) {
    return m2 * (s1.cross(s0).normalize());
  }

  /** Class representing XDS coordinate system */
  class XdsCoordinateSystem {

  public:

    /**
     * Initialise coordinate system. s0 should be the same length as s1, m2
     * should be a unit vector. These quantities are not checked because this
     * class will be created for each reflection and we want to maximize
     * performance.
     * @param s0 The incident beam vector
     * @param s1 The diffracted beam vector
     * @param m2 The rotation axis
     * @param phi The rotation angle
     */
    XdsCoordinateSystem(vec3 <double> s0,
                        vec3 <double> s1,
                        vec3 <double> m2,
                        double phi)
      : e1_(s1.cross(s0).normalize()),
        e2_(s1.cross(e1_).normalize()),
        e3_((s1 + s0).normalize()),
        zeta_(m2 * e1_) {}

  public:

    /** Get the e1 axis vector */
    vec3 <double> get_e1_axis() {
      return e1_;
    }

    /** Get the e2 axis vector */
    vec3 <double> get_e2_axis() {
      return e2_;
    }

    /** Get the e3 axis vector */
    vec3 <double> get_e3_axis() {
      return e3_;
    }

    /** Get the lorentz correction factor (zeta) */
    double get_zeta() {
      return zeta_;
    }

  private:

    vec3 <double> e1_;
    vec3 <double> e2_;
    vec3 <double> e3_;
    double zeta_;
  };

}} // namespace = dials::algorithms

#endif // DIALS_ALGORITHMS_INTEGRATION_XDS_COORDINATE_SYSTEM_H
