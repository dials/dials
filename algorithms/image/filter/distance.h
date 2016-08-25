/*
 * distance.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_IMAGE_FILTER_DISTANCE_H
#define DIALS_ALGORITHMS_IMAGE_FILTER_DISTANCE_H

#include <dials/array_family/scitbx_shared_and_versa.h>

namespace dials { namespace algorithms {

  /**
   * Compute the manhatten distance transform of a binary image
   * @param data The image
   * @return The distance transform
   */
  inline
  af::versa< int, af::c_grid<2> > manhatten_distance(
      const af::const_ref< bool, af::c_grid<2> > &data) {

    // Initialise stuff
    std::size_t height = data.accessor()[0];
    std::size_t width = data.accessor()[1];
    int max_distance = height + width;
    af::versa< int, af::c_grid<2> > result(data.accessor(), max_distance);

    // Go south and west
    for (std::size_t j = 0; j < height; ++j) {
      for (std::size_t i = 1; i < width; ++i) {
        int N = (j > 0) ? result(j-1,i) : max_distance;
        int E = (i > 0) ? result(j,i-1) : max_distance;
        if (data(j,i)) {
          result(j,i) = 0;
        } else {
          result(j,i) = 1 + std::min(N, E);
        }
      }
    }

    // Go north and east
    for (std::size_t j = height; j > 0; --j) {
      for (std::size_t i = width; i > 0; --i) {
        int S = (j < height) ? result(j,i-1) : max_distance;
        int W = (i < width)  ? result(j-1,i) : max_distance;
        int other = std::min(S, W);
        if (other < result(j-1,i-1)) {
          result(j-1,i-1) = 1 + other;
        }
      }
    }

    return result;
  }

}}

#endif // DIALS_ALGORITHMS_IMAGE_FILTER_DISTANCE_H
