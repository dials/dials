/*
 * rtmats.h
 *
 *  Copyright (C) (2016) STFC Rutherford Appleton Laboratory, UK.
 *
 *  Author: David Waterman.
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_REFINEMENT_RTMATS_H
#define DIALS_REFINEMENT_RTMATS_H

#ifndef DEG2RAD
#define DEG2RAD(x) ((x)*0.01745329251994329575)
#endif

#include <cmath>
#include <scitbx/mat3.h>
#include <scitbx/vec3.h>
#include <dials/error.h>

namespace dials { namespace refinement {

  using scitbx::mat3;
  using scitbx::vec3;

  /**
   * Calculate the first derivative of a rotation matrix R with respect to the
   * angle of rotation, given the axis and angle.
   *
   * Based on the FORTRAN subroutine RTMATS by David Thomas found in Mosflm.
   * Here the rotation is taken to be in a right-handed sense around the axis
   * whereas RTMATS uses a left-handed rotation.
   */
  mat3<double> dR_from_axis_and_angle(const vec3<double> &axis,
                                      double angle,
                                      bool deg = false) {
    if (deg) angle = DEG2RAD(angle);
    vec3<double> axis_ = axis.normalize();
    double ca = cos(angle);
    double sa = sin(angle);
    return mat3<double>(sa * axis_[0] * axis_[0] - sa,
                        sa * axis_[0] * axis_[1] - ca * axis_[2],
                        sa * axis_[0] * axis_[2] + ca * axis_[1],
                        sa * axis_[1] * axis_[0] + ca * axis_[2],
                        sa * axis_[1] * axis_[1] - sa,
                        sa * axis_[1] * axis_[2] - ca * axis_[0],
                        sa * axis_[2] * axis_[0] - ca * axis_[1],
                        sa * axis_[2] * axis_[1] + ca * axis_[0],
                        sa * axis_[2] * axis_[2] - sa);
  }

}}  // namespace dials::refinement

#endif  // DIALS_REFINEMENT_RTMATS_H
