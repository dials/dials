/*
 * unimodal.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_IMAGE_THRESHOLD_UNIMODAL_H
#define DIALS_ALGORITHMS_IMAGE_THRESHOLD_UNIMODAL_H

#include <cmath>
#include <iostream>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/ref_reductions.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using scitbx::af::int2;
  using scitbx::af::flex_int;
  using scitbx::af::flex_double;
  using scitbx::af::const_ref;
  using scitbx::af::max_index;
  using scitbx::af::max;
  using scitbx::af::min;

  /**
   * Calculate the maximum_deviation unimodal histogram threshold.
   * @param histo The histogram
   * @returns The threshold value
   */
  std::size_t maximum_deviation(const flex_double &histo) {

    // Get x, y at peak and at end of tail.
    std::size_t i0 = max_index(const_ref<double>(histo.begin(), histo.size()));
    std::size_t i1 = histo.size() - 1;
    double x0 = (double)i0 + 0.5, y0 = histo[i0];
    double x1 = (double)i1 + 0.5, y1 = histo[i1];

    // Calculate the line parameters
    double m = (y1 - y0) / (x1 - x0);
    double c = y0 - m * x0;

    // Find the maximum deviation from the line
    std::size_t imax = i0;
    double dmax = 0;
    for (std::size_t i = i0 + 1; i <= i1; ++i) {
      double x = i + 0.5;
      double y = histo[i];
      double d = std::abs(m * x - y + c) / std::sqrt(m * m + 1);
      if (d > dmax) {
        dmax = d;
        imax = i;
      }
    }

    // Return the maximum index
    return imax;
  }

  /**
   * Calculate the probability distribution of an image histogram
   * @param image The image to process
   * @param range The range of values to consider
   * @returns The probability distribution of values
   */
  flex_double probability_distribution(const flex_int &image, int2 range) {

    // Get the histogram range
    int minh = range[0];
    int maxi = max(const_ref<int>(image.begin(), image.size()));
    int maxh = min(int2(maxi, range[1]).const_ref());

    // Histogram the image
    flex_double p(maxh - minh + 1);
    std::size_t count = 0;
    for (std::size_t i = 0; i < image.size(); ++i) {
      if (minh <= image[i] && image[i] <= maxh) {
        p[image[i]-minh] += 1;
        count++;
      }
    }

    // Ensure more than 0 counts
    DIALS_ASSERT(count > 0);

    // Make into probability distribution
    for (std::size_t i = 0; i < p.size(); ++i) {
      p[i] /= count;
    }

    // Return the probability distribution
    return p;
  }

}} // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_IMAGE_THRESHOLD_UNIMODAL_H */
