/*
 * index_of_dispersion_filter.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_IMAGE_FILTER_INDEX_OF_DISPERSION_FILTER_H
#define DIALS_ALGORITHMS_IMAGE_FILTER_INDEX_OF_DISPERSION_FILTER_H

#include <algorithm>
#include <cmath>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/ref_reductions.h>
#include <dials/error.h>
#include "mean_and_variance.h"

namespace dials { namespace algorithms {

  using scitbx::af::int2;

  /**
   * Calculate the index of dispersion filtered image. The filtered image is
   * created by calculating the index of dispersion factor (variation index
   * var/mean) for each pixel.
   */
  template <typename FloatType = double>
  class IndexOfDispersionFilter {
  public:
    typedef FloatType value_type;

    /**
     * Calculate the index of dispersion filtered image. The filtered image is
     * created by calculating the index of dispersion factor (variation index
     * var/mean) for each pixel.
     * @param image The image to filter
     * @param size Size of the filter kernel (2 * size + 1)
     */
    IndexOfDispersionFilter(const af::const_ref<FloatType, af::c_grid<2> > &image,
                            int2 size) {
      // Get the mean and variance maps
      MeanAndVarianceFilter<FloatType> filter(image, size);
      mean_ = filter.mean();
      var_ = filter.sample_variance();

      // Calculate the filtered image
      index_of_dispersion_ = af::versa<FloatType, af::c_grid<2> >(var_.accessor(), 0);
      for (std::size_t i = 0; i < var_.size(); ++i) {
        if (mean_[i] > 0) {
          index_of_dispersion_[i] = var_[i] / mean_[i];
        } else {
          index_of_dispersion_[i] = 1.0;
        }
      }
    }

    /**
     * @returns The filtered image
     */
    af::versa<FloatType, af::c_grid<2> > index_of_dispersion() const {
      return index_of_dispersion_;
    }

    /**
     * @returns The mean filtered image
     */
    af::versa<FloatType, af::c_grid<2> > mean() const {
      return mean_;
    }

    /**
     * @returns The sample variance filtered image
     */
    af::versa<FloatType, af::c_grid<2> > sample_variance() const {
      return var_;
    }

  private:
    af::versa<FloatType, af::c_grid<2> > index_of_dispersion_;
    af::versa<FloatType, af::c_grid<2> > mean_;
    af::versa<FloatType, af::c_grid<2> > var_;
  };

  /**
   * Calculate the masked index of dispersion filtered image. The filtered image is
   * created by calculating the index of dispersion factor (variation index var/mean)
   * for each pixel. The mask is updated for those pixels where the filtered image is
   * invalid.
   */
  template <typename FloatType = double>
  class IndexOfDispersionFilterMasked {
  public:
    typedef FloatType value_type;

    /**
     * Calculate the masked index of dispersion filtered image. The filtered image is
     * created by calculating the index of dispersion factor (variation index var/mean)
     * for each pixel. The mask is updated for those pixels where the filtered image is
     * invalid.
     * @param image The image to filter
     * @param mask The mask to use
     * @param size Size of the filter kernel (2 * size + 1)
     * @param min_count The minimum counts under the filter to include the pixel
     */
    IndexOfDispersionFilterMasked(const af::const_ref<FloatType, af::c_grid<2> > &image,
                                  const af::const_ref<int, af::c_grid<2> > &mask,
                                  int2 size,
                                  int min_count) {
      // Get the mean and variance maps
      MeanAndVarianceFilterMasked<FloatType> filter(image, mask, size, min_count);
      mean_ = filter.mean();
      var_ = filter.sample_variance();
      mask_ = filter.mask();
      count_ = filter.count();

      // Calculate the filtered image.
      index_of_dispersion_ = af::versa<FloatType, af::c_grid<2> >(var_.accessor(), 0);
      for (std::size_t i = 0; i < var_.size(); ++i) {
        if (mask_[i] && mean_[i] > 0) {
          index_of_dispersion_[i] = var_[i] / mean_[i];
        } else {
          index_of_dispersion_[i] = 1.0;
          mask_[i] = 0;
        }
      }
    }

    /**
     * @returns The filter mask
     */
    af::versa<int, af::c_grid<2> > mask() const {
      return mask_;
    }

    /**
     * @returns The filter counts
     */
    af::versa<int, af::c_grid<2> > count() const {
      return count_;
    }

    /**
     * @returns The filtered image
     */
    af::versa<FloatType, af::c_grid<2> > index_of_dispersion() const {
      return index_of_dispersion_;
    }

    /**
     * @returns The mean filtered image
     */
    af::versa<FloatType, af::c_grid<2> > mean() const {
      return mean_;
    }

    /**
     * @returns The sample variance filtered image
     */
    af::versa<FloatType, af::c_grid<2> > sample_variance() const {
      return var_;
    }

  private:
    af::versa<int, af::c_grid<2> > count_;
    af::versa<int, af::c_grid<2> > mask_;
    af::versa<FloatType, af::c_grid<2> > index_of_dispersion_;
    af::versa<FloatType, af::c_grid<2> > mean_;
    af::versa<FloatType, af::c_grid<2> > var_;
  };

}}  // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_IMAGE_FILTER_INDEX_OF_DISPERSION_FILTER_H */
