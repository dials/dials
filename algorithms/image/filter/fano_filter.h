/*
 * fano_filter.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_IMAGE_FILTER_FANO_FILTER_H
#define DIALS_ALGORITHMS_IMAGE_FILTER_FANO_FILTER_H

#include <algorithm>
#include <cmath>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/ref_reductions.h>
#include <dials/error.h>
#include "mean_and_variance.h"

namespace dials { namespace algorithms {

  using scitbx::af::int2;
  using scitbx::af::flex;
  using scitbx::af::flex_int;
  using scitbx::af::flex_bool;
  using scitbx::af::flex_double;


  /**
   * Calculate the fano filtered image. The filtered image is created by
   * calculating the fano factor (variation index var/mean) for each pixel.
   */
  class FanoFilter {
  public:


    /**
     * Calculate the fano filtered image. The filtered image is created by
     * calculating the fano factor (variation index var/mean) for each pixel.
     * @param image The image to filter
     * @param size Size of the filter kernel (2 * size + 1)
     */
    FanoFilter(const flex_double &image, int2 size) {

      // Get the mean and variance maps
      MeanAndVarianceFilter filter(image, size);
      mean_ = filter.mean();
      var_  = filter.sample_variance();

      // Calculate the filtered image
      fano_ = flex_double(var_.accessor(), 0);
      for (std::size_t i = 0; i < var_.size(); ++i) {
        if (mean_[i] > 0) {
          fano_[i] = var_[i] / mean_[i];
        } else {
          fano_[i] = 1.0;
        }
      }
    }

    /**
     * @returns The filter mask
     */
    flex_int mask() const {
      return mask_;
    }

    /**
     * @returns The filtered image
     */
    flex_double fano() const {
      return fano_;
    }

    /**
     * @returns The mean filtered image
     */
    flex_double mean() const {
      return mean_;
    }

    /**
     * @returns The sample variance filtered image
     */
    flex_double sample_variance() const {
      return var_;
    }

  private:
    flex_int mask_;
    flex_double fano_;
    flex_double mean_;
    flex_double var_;
  };


  /**
   * Calculate the masked fano filtered image. The filtered image is created by
   * calculating the fano factor (variation index var/mean) for each pixel.
   * The mask is updated for those pixels where the filtered image is invalid.
   */
  class FanoFilterMasked {
  public:


    /**
     * Calculate the masked fano filtered image. The filtered image is created
     * by calculating the fano factor (variation index var/mean) for each pixel.
     * The mask is updated for those pixels where the filtered image is invalid.
     * @param image The image to filter
     * @param mask The mask to use
     * @param size Size of the filter kernel (2 * size + 1)
     * @param min_count The minimum counts under the filter to include the pixel
     */
    FanoFilterMasked(const flex_double &image, const flex_int &mask,
        int2 size, int min_count) {

      // Get the mean and variance maps
      MeanAndVarianceFilterMasked filter(image, mask, size, min_count);
      mean_  = filter.mean();
      var_   = filter.sample_variance();
      mask_  = filter.mask();
      count_ = filter.count();

      // Calculate the filtered image.
      fano_ = flex_double(var_.accessor(), 0);
      for (std::size_t i = 0; i < var_.size(); ++i) {
        if (mask_[i] && mean_[i] > 0) {
          fano_[i] = var_[i] / mean_[i];
        } else {
          fano_[i] = 1.0;
          mask_[i] = 0;
        }
      }
    }

    /**
     * @returns The filter mask
     */
    flex_int mask() const {
      return mask_;
    }

    /**
     * @returns The filter counts
     */
    flex_int count() const {
      return count_;
    }

    /**
     * @returns The filtered image
     */
    flex_double fano() const {
      return fano_;
    }

    /**
     * @returns The mean filtered image
     */
    flex_double mean() const {
      return mean_;
    }

    /**
     * @returns The sample variance filtered image
     */
    flex_double sample_variance() const {
      return var_;
    }

  private:
    flex_int count_;
    flex_int mask_;
    flex_double fano_;
    flex_double mean_;
    flex_double var_;
  };

}} // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_IMAGE_FILTER_FANO_FILTER_H */
