/*
 * simple.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_IMAGE_FILL_HOLES_SIMPLE_H
#define DIALS_ALGORITHMS_IMAGE_FILL_HOLES_SIMPLE_H

#include <vector>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/algorithms/image/filter/distance.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  namespace detail {

    /**
     * Struct to provide a pixel node
     */
    struct SimpleFillNode {
      std::size_t j;
      std::size_t i;
      int d;
      SimpleFillNode(std::size_t j_, std::size_t i_, int d_) : j(j_), i(i_), d(d_) {}
    };

    /**
     * Struct to provide pixel node comparison
     */
    struct CompareSimpleFillNode {
      bool operator()(const SimpleFillNode &a, const SimpleFillNode &b) const {
        return a.d < b.d;
      }
    };
  }  // namespace detail

  /**
   * A simple function to fill holes in images
   * @param data The data array
   * @param mask The mask array
   * @returns The filled image
   */
  inline af::versa<double, af::c_grid<2> > simple_fill(
    const af::const_ref<double, af::c_grid<2> > &data,
    const af::const_ref<bool, af::c_grid<2> > &mask) {
    // Check the input
    DIALS_ASSERT(data.accessor().all_eq(mask.accessor()));
    std::size_t height = data.accessor()[0];
    std::size_t width = data.accessor()[1];

    // Compute the manhattan distance transform of the mask
    af::versa<int, af::c_grid<2> > distance(mask.accessor());
    manhattan_distance(mask, true, distance.ref());

    // Get a list of the pixels to fill
    std::vector<detail::SimpleFillNode> pixels;
    for (std::size_t j = 0; j < height; ++j) {
      for (std::size_t i = 0; i < width; ++i) {
        if (distance(j, i) > 0) {
          pixels.push_back(detail::SimpleFillNode(j, i, distance(j, i)));
        }
      }
    }

    // Sort in order of distance
    std::sort(pixels.begin(), pixels.end(), detail::CompareSimpleFillNode());

    // Fill in pixels
    af::versa<double, af::c_grid<2> > result(data.accessor());
    std::copy(data.begin(), data.end(), result.begin());
    for (std::size_t k = 0; k < pixels.size(); ++k) {
      std::size_t j = pixels[k].j;
      std::size_t i = pixels[k].i;
      double sum = 0.0;
      std::size_t num = 0;
      if (j > 0 && distance(j - 1, i) == 0) {
        sum += result(j - 1, i);
        num += 1;
      }
      if (i > 0 && distance(j, i - 1) == 0) {
        sum += result(j, i - 1);
        num += 1;
      }
      if (j < height - 1 && distance(j + 1, i) == 0) {
        sum += result(j + 1, i);
        num += 1;
      }
      if (i < width - 1 && distance(j, i + 1) == 0) {
        sum += result(j, i + 1);
        num += 1;
      }
      DIALS_ASSERT(num > 0);
      result(j, i) = sum / (double)num;
      distance(j, i) = 0;
    }

    // Ensure every pixel has been used
    DIALS_ASSERT(distance.all_eq(0));

    // Return the filled pixels
    return result;
  }

  /**
   * A simple function to fill holes in images
   * @param data The data array
   * @param mask The mask array
   * @returns The filled image
   */
  inline af::versa<double, af::c_grid<2> > diffusion_fill(
    const af::const_ref<double, af::c_grid<2> > &data,
    const af::const_ref<bool, af::c_grid<2> > &mask,
    std::size_t niter) {
    // Check input
    DIALS_ASSERT(niter > 0);
    DIALS_ASSERT(data.accessor().all_eq(mask.accessor()));
    std::size_t height = data.accessor()[0];
    std::size_t width = data.accessor()[1];

    // Copy initial values
    af::versa<double, af::c_grid<2> > result(data.accessor());
    for (std::size_t i = 0; i < data.size(); ++i) {
      result[i] = data[i];
    }

    // Fill missing values
    for (std::size_t iter = 0; iter < niter; ++iter) {
      for (std::size_t j = 0; j < height; ++j) {
        for (std::size_t i = 0; i < width; ++i) {
          if (!mask(j, i)) {
            double sum = 0;
            int cnt = 0;
            if (i > 0) {
              sum += result(j, i - 1);
              cnt += 1;
            }
            if (j > 0) {
              sum += result(j - 1, i);
              cnt += 1;
            }
            if (i < width - 1) {
              sum += result(j, i + 1);
              cnt += 1;
            }
            if (j < height - 1) {
              sum += result(j + 1, i);
              cnt += 1;
            }
            result(j, i) = sum / (double)cnt;
          }
        }
      }
    }

    return result;
  }

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_IMAGE_FILL_HOLES_SIMPLE_H
