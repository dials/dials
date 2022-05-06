/*
 * median.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_IMAGE_FILTER_MEDIAN_H
#define DIALS_ALGORITHMS_IMAGE_FILTER_MEDIAN_H

#include <algorithm>
#include <scitbx/array_family/tiny_types.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using scitbx::af::int2;

  /**
   * Apply a median filter to an image
   * @param image The image to filter
   * @param size The size of the filter kernel
   * @returns The filtered image
   */
  template <typename T>
  af::versa<T, af::c_grid<2> > median_filter(
    const af::const_ref<T, af::c_grid<2> > &image,
    int2 size) {
    // Check the input is valid
    DIALS_ASSERT(size.all_ge(0));
    DIALS_ASSERT(image.accessor().all_gt(0));

    // The array for output
    af::versa<T, af::c_grid<2> > median(image.accessor(), T(0));

    // Create the array to sort to get the median
    std::size_t ysize = image.accessor()[0];
    std::size_t xsize = image.accessor()[1];
    af::shared<T> pixels_array((2 * size[0] + 1) * (2 * size[1] + 1), T(0));
    af::ref<T> pixels = pixels_array.ref();

    // For each pixel obtain the median value
    for (std::size_t j = 0; j < ysize; ++j) {
      for (std::size_t i = 0; i < xsize; ++i) {
        std::size_t npix = 0;
        for (int jj = (int)j - size[0]; jj <= (int)j + size[0]; ++jj) {
          for (int ii = (int)i - size[1]; ii <= (int)i + size[1]; ++ii) {
            if (jj >= 0 && ii >= 0 && jj < ysize && ii < xsize) {
              pixels[npix++] = image(jj, ii);
            }
          }
        }
        size_t n = npix / 2;
        std::nth_element(pixels.begin(), pixels.begin() + n, pixels.begin() + npix);
        median(j, i) = pixels[n];
      }
    }

    // Return the median filtered image
    return median;
  }

  /**
   * Apply a median filter to an image with a mask
   * @param image The image to filter
   * @param mask The image mask
   * @param size The size of the filter kernel
   * @param periodic Wrap the filter
   * @returns The filtered image
   */
  template <typename T>
  af::versa<T, af::c_grid<2> > median_filter_masked(
    const af::const_ref<T, af::c_grid<2> > &image,
    const af::const_ref<bool, af::c_grid<2> > &mask,
    int2 size,
    bool periodic) {
    // Check the input is valid
    DIALS_ASSERT(size.all_ge(0));
    DIALS_ASSERT(image.accessor().all_gt(0));
    DIALS_ASSERT(image.accessor().all_eq(mask.accessor()));

    // The array for output
    af::versa<T, af::c_grid<2> > median(image.accessor(), T(0));

    // Create the array to sort to get the median
    std::size_t ysize = image.accessor()[0];
    std::size_t xsize = image.accessor()[1];
    af::shared<T> pixels_array((2 * size[0] + 1) * (2 * size[1] + 1), T(0));
    af::ref<T> pixels = pixels_array.ref();

    // For each pixel obtain the median value
    for (int j = 0; j < (int)ysize; ++j) {
      for (int i = 0; i < (int)xsize; ++i) {
        std::size_t npix = 0;
        for (int jj = j - size[0]; jj <= j + size[0]; ++jj) {
          for (int ii = i - size[1]; ii <= i + size[1]; ++ii) {
            if (periodic) {
              std::size_t jjj = (jj % ysize + ysize) % ysize;
              std::size_t iii = (ii % xsize + xsize) % xsize;
              DIALS_ASSERT(jjj >= 0 && iii >= 0 && jjj < ysize && iii < xsize);
              if (mask(jjj, iii)) {
                DIALS_ASSERT(npix < pixels.size());
                pixels[npix++] = image(jjj, iii);
              }
            } else {
              if (jj >= 0 && ii >= 0 && jj < ysize && ii < xsize) {
                if (mask(jj, ii)) {
                  pixels[npix++] = image(jj, ii);
                }
              }
            }
          }
        }
        DIALS_ASSERT(npix <= pixels.size());
        if (npix > 0) {
          size_t n = npix / 2;
          std::nth_element(pixels.begin(), pixels.begin() + n, pixels.begin() + npix);
          median(j, i) = pixels[n];
        }
      }
    }

    // Return the median filtered image
    return median;
  }

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_IMAGE_FILTER_MEDIAN_H
