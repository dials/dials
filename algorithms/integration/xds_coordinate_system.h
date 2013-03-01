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

  /** Class representing XDS coordinate system */
  class XdsCoordinateSystem {

  public:

    /**
     * Initialise coordinate system
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
