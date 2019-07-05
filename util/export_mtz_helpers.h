/*
 * FIXME add a header
 */

#ifndef DIALS_UTIL_EXPORT_MTZ_HELPERS_H
#define DIALS_UTIL_EXPORT_MTZ_HELPERS_H

#include <cmath>
#include <scitbx/vec3.h>
#include <scitbx/vec2.h>
#include <scitbx/constants.h>
#include <scitbx/matrix/multiply.h>
#include <scitbx/math/r3_rotation.h>
#include <scitbx/array_family/tiny_types.h>
#include <cctbx/uctbx.h>

namespace dials { namespace util {

  using cctbx::uctbx::unit_cell;
  using scitbx::mat3;

  mat3<double> dials_u_to_mosflm(const mat3<double> dials_U, unit_cell uc) {
    scitbx::af::double6 p = uc.parameters();
    scitbx::af::double6 rp = uc.reciprocal_parameters();
    scitbx::mat3<double> dials_B = uc.fractionalization_matrix().transpose();
    scitbx::mat3<double> dials_UB = dials_U * dials_B;

    double d2r = atan(1.0) / 45.0;

    scitbx::mat3<double> mosflm_B(rp[0],
                                  rp[1] * cos(d2r * rp[5]),
                                  rp[2] * cos(d2r * rp[4]),
                                  0,
                                  rp[1] * sin(d2r * rp[5]),
                                  -rp[2] * sin(d2r * rp[4]) * cos(d2r * p[3]),
                                  0,
                                  0,
                                  1.0 / p[2]);

    scitbx::mat3<double> mosflm_U = dials_UB * mosflm_B.inverse();

    return mosflm_U;
  }
}}  // namespace dials::util

#endif
