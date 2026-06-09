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
#include <iotbx/mtz/object.h>
#include <scitbx/array_family/misc_functions.h>

namespace dials { namespace util {

  using cctbx::uctbx::unit_cell;
  using iotbx::mtz::object;
  using scitbx::mat3;
  using scitbx::af::min_index;

  mat3<double> ub_to_mosflm_u(const mat3<double> UB, unit_cell uc) {
    scitbx::af::double6 p = uc.parameters();
    scitbx::af::double6 rp = uc.reciprocal_parameters();

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

    scitbx::mat3<double> mosflm_U = UB * mosflm_B.inverse();

    return mosflm_U;
  }

}}  // namespace dials::util

#endif
