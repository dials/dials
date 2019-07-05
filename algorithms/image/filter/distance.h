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
#include <dials/error.h>

namespace dials { namespace algorithms {

  /**
   * Compute the manhattan distance transform of a binary image
   * @param data The image
   * @return The distance transform
   */
  template <typename InputType, typename OutputType>
  void manhattan_distance(const af::const_ref<InputType, af::c_grid<2> > &src,
                          InputType value,
                          af::ref<OutputType, af::c_grid<2> > dst) {
    // Initialise stuff
    std::size_t height = src.accessor()[0];
    std::size_t width = src.accessor()[1];
    OutputType max_distance = height + width;
    DIALS_ASSERT(src.accessor().all_eq(dst.accessor()));

    // Go north and east
    for (std::size_t j = 0; j < height; ++j) {
      for (std::size_t i = 1; i < width; ++i) {
        OutputType N = (j > 0) ? dst(j - 1, i) : max_distance;
        OutputType E = (i > 0) ? dst(j, i - 1) : max_distance;
        if (src(j, i) == value) {
          dst(j, i) = 0;
        } else {
          dst(j, i) = 1 + std::min(N, E);
        }
      }
    }

    // Go south and west
    for (std::size_t j = height; j > 0; --j) {
      for (std::size_t i = width; i > 0; --i) {
        OutputType S = (j < height) ? dst(j, i - 1) : max_distance;
        OutputType W = (i < width) ? dst(j - 1, i) : max_distance;
        OutputType other = std::min(S, W);
        if (other < dst(j - 1, i - 1)) {
          dst(j - 1, i - 1) = 1 + other;
        }
      }
    }
  }

  /**
   * Compute the chebyshev distance to the given value
   * @param src: The src array
   * @param value: The value to compute the distance to
   * @param dst: The destination array
   */
  template <typename InputType, typename OutputType>
  void chebyshev_distance(const af::const_ref<InputType, af::c_grid<2> > &src,
                          InputType value,
                          af::ref<OutputType, af::c_grid<2> > dst) {
    // Initialise stuff
    std::size_t height = src.accessor()[0];
    std::size_t width = src.accessor()[1];
    OutputType max_distance = height + width;
    DIALS_ASSERT(src.accessor().all_eq(dst.accessor()));

    // Go north and east
    for (std::size_t j = 0; j < height; ++j) {
      for (std::size_t i = 1; i < width; ++i) {
        OutputType N = (j > 0) ? dst(j - 1, i) : max_distance;
        OutputType E = (i > 0) ? dst(j, i - 1) : max_distance;
        OutputType NE = (j > 0 && i > 0) ? dst(j - 1, i - 1) : max_distance;
        OutputType NW = (j > 0 && i < width - 1) ? dst(j - 1, i + 1) : max_distance;
        if (src(j, i) == value) {
          dst(j, i) = 0;
        } else {
          dst(j, i) = 1 + std::min(std::min(N, E), std::min(NE, NW));
        }
      }
    }

    // Go south and west
    for (std::size_t j = height; j > 0; --j) {
      for (std::size_t i = width; i > 0; --i) {
        OutputType S = (j < height) ? dst(j, i - 1) : max_distance;
        OutputType W = (i < width) ? dst(j - 1, i) : max_distance;
        OutputType SE = (j < height && i > 1) ? dst(j, i - 2) : max_distance;
        OutputType SW = (j < height && i < width) ? dst(j, i) : max_distance;
        OutputType other = std::min(std::min(S, W), std::min(SE, SW));
        if (other < dst(j - 1, i - 1)) {
          dst(j - 1, i - 1) = 1 + other;
        }
      }
    }
  }
}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_IMAGE_FILTER_DISTANCE_H
