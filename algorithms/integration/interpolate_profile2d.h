

#ifndef DIALS_ALGORITHMS_INTEGRATION_INERPOLATE_PROFILE_2D_H
#define DIALS_ALGORITHMS_INTEGRATION_INERPOLATE_PROFILE_2D_H

#include <scitbx/vec3.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>


namespace dials { namespace algorithms {

  using scitbx::vec3;

  af::versa< double, af::c_grid<2> > interpolate_profile2d(
      double x, double y,
      const af::const_ref< double, af::c_grid<2> > &data) {
    int ox = 0, oy = 0;
    DIALS_ASSERT(x >= -0.5 && x <= 0.5);
    DIALS_ASSERT(y >= -0.5 && y <= 0.5);
    if (x < 0) {
      ox = -1;
      x += 1;
    }
    if (y < 0) {
      oy = -1;
      y += 1;
    }
    std::size_t h = data.accessor()[0];
    std::size_t w = data.accessor()[1];
    af::versa< double, af::c_grid<2> > result(data.accessor());
    for (std::size_t j = 0; j < h; ++j) {
      for (std::size_t i = 0; i < w; ++i) {
        int j0 = j + oy;
        int i0 = i + ox;
        int j1 = j0 + 1;
        int i1 = i0 + 1;
        if (j0 < 0) j0 = 0;
        if (i0 < 0) i0 = 0;
        if (j1 >= h) j1 = h-1;
        if (i1 >= w) i1 = w-1;
        double f00 = data(j0,i0);
        double f10 = data(j1,i0);
        double f01 = data(j0,i1);
        double f11 = data(j1,i1);
        double I = f00*(1-x)*(1-y)+f10*y*(1-x)+f01*x*(1-y)+f11*x*y;
        result(j,i) = I;
      }
    }
    return result;
  }

}}

#endif // DIALS_ALGORITHMS_INTEGRATION_INERPOLATE_PROFILE_2D_H
