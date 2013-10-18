/*
 * ideal_profile.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_REFLECTION_BASIS_IDEAL_PROFILE_H
#define DIALS_ALGORITHMS_REFLECTION_BASIS_IDEAL_PROFILE_H

#include <cmath>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <scitbx/array_family/ref_reductions.h>

namespace dials { namespace algorithms { namespace reflection_basis {

  /**
   * Evalaute a gaussian at a point
   * @param x The coordinate
   * @param x0 The centre of the distribution
   * @param s The standard deviation of the distribution
   * @returns The value of the gaussian at that point
   */
  inline
  double evaluate_gaussian(double x, double x0, double s) {
    double v = (x - x0) / s;
    return std::exp(-(v * v) / 2.0);
  }

  /**
   * Generate an ideal profile in the reflection frame
   * @param size The size of the grid (2 * size + 1)
   * @param nsig The number of standard deviations
   * @returns The profile
   */
  inline
  af::versa< double, af::c_grid<3> >
  ideal_profile(std::size_t size, std::size_t nsig) {
    double centre = size;
    size = 2 * size + 1;
    double sig = centre / nsig;

    af::c_grid<3> accessor(size, size, size);
    af::versa< double, af::c_grid<3> > profile(accessor, 0.0);
    for (std::size_t k = 0; k < size; ++k) {
      for (std::size_t j = 0; j < size; ++j) {
        for (std::size_t i = 0; i < size; ++i) {
          double gx = evaluate_gaussian(i, centre, sig);
          double gy = evaluate_gaussian(j, centre, sig);
          double gz = evaluate_gaussian(k, centre, sig);
          profile(k, j, i) = gx * gy * gz;
        }
      }
    }
    double tot = scitbx::af::sum(profile.const_ref());
    DIALS_ASSERT(tot > 0);
    for (std::size_t i = 0; i < profile.size(); ++i) {
      profile[i] /= tot;
    }
    return profile;
  }

}}} // namespace dials::algorithms::reflection_basis

#endif /* DIALS_ALGORITHMS_REFLECTION_BASIS_IDEAL_PROFILE_H */
