/*
 * convolve.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_IMAGE_FILTER_CONVOLVE_H
#define DIALS_ALGORITHMS_IMAGE_FILTER_CONVOLVE_H

#include <omptbx/omp_or_stubs.h>

#include <algorithm>
#include <cmath>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/ref_reductions.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using scitbx::af::int2;

  /**
   * Perform a simple convolution between an image and kernel
   * @param image The image to filter
   * @param kernel The kernel to convolve with
   * @returns The convolved image
   */
  template <typename FloatType>
  af::versa<FloatType, af::c_grid<2> > convolve(
    const af::const_ref<FloatType, af::c_grid<2> > &image,
    const af::const_ref<FloatType, af::c_grid<2> > &kernel) {
    // Only allow odd-sized kernel sizes
    DIALS_ASSERT(kernel.accessor()[0] & 1);
    DIALS_ASSERT(kernel.accessor()[1] & 1);

    // The image sizes and mid-point
    int2 isz = image.accessor();
    int2 ksz = kernel.accessor();
    int2 mid(ksz[0] / 2, ksz[1] / 2);

    // Create the output
    af::versa<FloatType, af::c_grid<2> > result(image.accessor(),
                                                af::init_functor_null<FloatType>());

    // Convolve the image with the kernel
    for (int j = 0; j < isz[0]; ++j) {
      for (int i = 0; i < isz[1]; ++i) {
        result(j, i) = 0.0;
        for (int jj = 0; jj < ksz[0]; ++jj) {
          for (int ii = 0; ii < ksz[1]; ++ii) {
            int jjj = j + jj - mid[0];
            int iii = i + ii - mid[1];
            if (jjj < 0)
              jjj = 0;
            else if (jjj >= isz[0])
              jjj = isz[0] - 1;
            if (iii < 0)
              iii = 0;
            else if (iii >= isz[1])
              iii = isz[1] - 1;
            result(j, i) += image(jjj, iii) * kernel(jj, ii);
          }
        }
      }
    }

    // Return the result
    return result;
  }

  /**
   * Perform a seperable row convolution between an image and kernel
   * @param image The image to filter
   * @param kernel The kernel to convolve with
   * @returns The convolved image
   */
  template <typename FloatType>
  af::versa<FloatType, af::c_grid<2> > convolve_row(
    const af::const_ref<FloatType, af::c_grid<2> > &image,
    const af::const_ref<FloatType> &kernel) {
    // Only allow odd-sized kernel sizes
    DIALS_ASSERT(kernel.size() & 1);

    // The image sizes and mid-point
    int2 isz = image.accessor();
    std::size_t ksz = kernel.size();
    std::size_t mid = ksz / 2;

    // Create the output
    af::versa<FloatType, af::c_grid<2> > result(image.accessor(),
                                                af::init_functor_null<FloatType>());

    // Convolve the image with the kernel
    for (int j = 0; j < isz[0]; ++j) {
      for (int i = 0; i < isz[1]; ++i) {
        result(j, i) = 0.0;
        for (int ii = 0; ii < ksz; ++ii) {
          int iii = i + ii - mid;
          if (iii < 0)
            iii = 0;
          else if (iii >= isz[1])
            iii = isz[1] - 1;
          result(j, i) += image(j, iii) * kernel[ii];
        }
      }
    }

    // Return the result
    return result;
  }

  /**
   * Perform a seperable column convolution between an image and kernel
   * @param image The image to filter
   * @param kernel The kernel to convolve with
   * @returns The convolved image
   */
  template <typename FloatType>
  af::versa<FloatType, af::c_grid<2> > convolve_col(
    const af::const_ref<FloatType, af::c_grid<2> > &image,
    const af::const_ref<FloatType> &kernel) {
    // Only allow odd-sized kernel sizes
    DIALS_ASSERT(kernel.size() & 1);

    // The image sizes and mid-point
    int2 isz = image.accessor();
    std::size_t ksz = kernel.size();
    std::size_t mid = ksz / 2;

    // Create the output
    af::versa<FloatType, af::c_grid<2> > result(image.accessor(),
                                                af::init_functor_null<FloatType>());

    // Convolve the image with the kernel
    for (int j = 0; j < isz[0]; ++j) {
      for (int i = 0; i < isz[1]; ++i) {
        result(j, i) = 0.0;
        for (int jj = 0; jj < ksz; ++jj) {
          int jjj = j + jj - mid;
          if (jjj < 0)
            jjj = 0;
          else if (jjj >= isz[0])
            jjj = isz[0] - 1;
          result(j, i) += image(jjj, i) * kernel[jj];
        }
      }
    }

    // Return the result
    return result;
  }

}}  // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_IMAGE_FILTER_CONVOLVE_H */
