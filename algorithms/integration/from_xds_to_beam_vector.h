/*
 * from_xds_to_beam_vector.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_INTEGRATION_FROM_XDS_TO_BEAM_VECTOR_H
#define DIALS_ALGORITHMS_INTEGRATION_FROM_XDS_TO_BEAM_VECTOR_H

#include <cmath>
#include <scitbx/vec3.h>
#include <dials/error.h>
#include "xds_coordinate_system.h"

namespace dials { namespace algorithms {

  using scitbx::vec3;

  /** A class to perform a transform from XDS coordinate frame to beam vector */
  class FromXdsToBeamVector {

  public:

    /**
     * Initialise the transform
     * @param xcs The XDS coordinate system
     * @param s1 The beam vector
     */
    FromXdsToBeamVector(XdsCoordinateSystem xcs,
                        vec3 <double> s1)
      : scaled_e1_(xcs.get_e1_axis() * s1.length()),
        scaled_e2_(xcs.get_e2_axis() * s1.length()),
        normalized_s1_(s1.normalize()),
        radius_(s1.length()) {}

    /**
     * Apply the transform to the xds point. The transform is done by
     * calculating the coordinate of point in the sphere tangent plane at the
     * location of the exit point of the beam vector, s1. The beam vector, s',
     * is then calculated by finding the intersection of line defined by the
     * plane normal (i.e. s1 vector) with origin at the calculated tangent plane
     * point point and the ewald sphere.
     *
     * @param c The XDS coordinate
     * @returns The beam vector
     */
    vec3 <double> operator()(vec3 <double> c) const {
      vec3 <double> p = c[0] * scaled_e1_ + c[1] * scaled_e2_;
      double b = radius_ * radius_ - p.length_sq();
      DIALS_ASSERT(b >= 0);
      double d = -(normalized_s1_ * p) + std::sqrt(b);
      return p + d * normalized_s1_;
    }

  private:

    vec3 <double> scaled_e1_;
    vec3 <double> scaled_e2_;
    vec3 <double> normalized_s1_;
    double radius_;
  };

}} // namespace = dials::algorithms

#endif // DIALS_ALGORITHMS_INTEGRATION_FROM_XDS_TO_BEAM_VECTOR_H
