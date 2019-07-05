/*
 * mean_and_variance.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_IMAGE_FILTER_MEAN_AND_VARIANCE_H
#define DIALS_ALGORITHMS_IMAGE_FILTER_MEAN_AND_VARIANCE_H

#include <algorithm>
#include <cmath>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/ref_reductions.h>
#include <dials/error.h>
#include "summed_area.h"

namespace dials { namespace algorithms {

  using scitbx::af::int2;

  /**
   * Calculate the mean box filtered image.
   * @param image The image to filter
   * @param size The size of the filter kernel (2 * size + 1)
   * @returns The filtered image
   */
  template <typename FloatType>
  af::versa<FloatType, af::c_grid<2> > mean_filter(
    const af::const_ref<FloatType, af::c_grid<2> > &image,
    int2 size) {
    // Check the input is valid
    DIALS_ASSERT(size.all_gt(0));
    DIALS_ASSERT(image.accessor().all_gt(0));

    // Get the summed area image
    af::versa<FloatType, af::c_grid<2> > mean = summed_area<FloatType>(image, size);
    FloatType inv_count = 1.0 / ((FloatType)(2 * size[0] + 1) * (2 * size[1] + 1));

    // Calculate the mean at each point
    for (std::size_t i = 0; i < mean.size(); ++i) {
      mean[i] *= inv_count;
    }
    return mean;
  }

  /**
   * Calculate the masked mean box filtered image. The minimum number of
   * counts specifies which pixels near a masked region will be used.
   *
   * If the minimum counts == 1 then all non-masked pixels will be used.
   * If the minimum counts == size of the kernel then only though pixels
   * where all surrounding pixels are on will be used etc.
   *
   * A min_count value <= 0 will be interpreted as the size of the kernel
   *
   * The mask is updated with those pixels which are not filtered.
   *
   * @param image The image to filter
   * @param mask The mask to use (0 = off, 1 == on)
   * @param size The size of the filter kernel (2 * size + 1)
   * @param min_count The minimum counts to use
   * @param ignore_masked Ignore and set mean to zero if masked
   * @returns The filtered image
   */
  template <typename FloatType>
  af::versa<FloatType, af::c_grid<2> > mean_filter_masked(
    const af::const_ref<FloatType, af::c_grid<2> > &image,
    af::ref<int, af::c_grid<2> > mask,
    int2 size,
    int min_count,
    bool ignore_masked = true) {
    // Check the input is valid
    DIALS_ASSERT(size.all_ge(0));
    DIALS_ASSERT(image.accessor().all_gt(0));
    DIALS_ASSERT(image.accessor().all_eq(mask.accessor()));

    // Ensure the min counts are valid
    if (min_count <= 0) {
      min_count = (2 * size[0] + 1) * (2 * size[1] + 1);
    } else {
      DIALS_ASSERT(min_count <= (2 * size[0] + 1) * (2 * size[1] + 1));
    }

    // Calculate the summed area under the mask
    af::versa<int, af::c_grid<2> > summed_mask = summed_area<int>(mask, size);

    // Ensure that all masked pixels are zero in the image and update the mask
    af::versa<FloatType, af::c_grid<2> > temp(image.accessor(),
                                              af::init_functor_null<FloatType>());
    for (std::size_t i = 0; i < image.size(); ++i) {
      temp[i] = image[i] * (mask[i] != 0);
      mask[i] *= (summed_mask[i] >= min_count);
    }

    // Calculate the summed area under the image
    af::versa<FloatType, af::c_grid<2> > summed_image =
      summed_area<FloatType>(temp.const_ref(), size);

    // Calculate the mean filtered image
    for (std::size_t i = 0; i < image.size(); ++i) {
      if ((!ignore_masked || mask[i]) && summed_mask[i] >= min_count) {
        summed_image[i] /= (FloatType)summed_mask[i];
      } else {
        summed_image[i] = 0.0;
      }
    }

    // Return the filtered image
    return summed_image;
  }

  /**
   * A class to calculate the mean and variance filtered images.
   */
  template <typename FloatType = double>
  class MeanAndVarianceFilter {
  public:
    typedef FloatType value_type;

    /**
     * Initialise the algorithm.
     * @params image The image to filter
     * @param size The size of the filter kernel (2 * size + 1)
     */
    MeanAndVarianceFilter(const af::const_ref<FloatType, af::c_grid<2> > &image,
                          int2 size) {
      // Check the input is valid
      DIALS_ASSERT(size.all_gt(0));
      DIALS_ASSERT(image.accessor().all_gt(0));

      // Calculate the squared image
      af::versa<FloatType, af::c_grid<2> > image_sq(image.accessor(),
                                                    af::init_functor_null<FloatType>());
      for (std::size_t i = 0; i < image.size(); ++i) {
        image_sq[i] = image[i] * image[i];
      }

      // Inverse of counts to avoid excessive division
      int count = (2 * size[0] + 1) * (2 * size[1] + 1);
      inv_count_ = 1.0 / count;
      inv_countm1_ = 1.0 / (count - 1);

      // Calculate the summed area under the image and image**2
      sum_ = summed_area<FloatType>(image, size);
      sq_sum_ = summed_area<FloatType>(image_sq.const_ref(), size);
    }

    /**
     * Calculate the mean filtered image
     * @returns The mean filtered image.
     */
    af::versa<FloatType, af::c_grid<2> > mean() const {
      af::versa<FloatType, af::c_grid<2> > m(sum_.accessor(),
                                             af::init_functor_null<FloatType>());
      for (std::size_t i = 0; i < sum_.size(); ++i) {
        m[i] = sum_[i] * inv_count_;
      }
      return m;
    }

    /**
     * Calculate the variance filtered image
     * @returns The variance filtered image
     */
    af::versa<FloatType, af::c_grid<2> > variance() const {
      af::versa<FloatType, af::c_grid<2> > v(sum_.accessor(),
                                             af::init_functor_null<FloatType>());
      for (std::size_t i = 0; i < sum_.size(); ++i) {
        v[i] = (sq_sum_[i] - (sum_[i] * sum_[i] * inv_count_)) * (inv_count_);
      }
      return v;
    }

    /**
     * Calculate the sample variance filtered image
     * @returns The sample variance filtered image
     */
    af::versa<FloatType, af::c_grid<2> > sample_variance() const {
      af::versa<FloatType, af::c_grid<2> > v(sum_.accessor(),
                                             af::init_functor_null<FloatType>());
      for (std::size_t i = 0; i < sum_.size(); ++i) {
        v[i] = (sq_sum_[i] - (sum_[i] * sum_[i] * inv_count_)) * (inv_countm1_);
      }
      return v;
    }

  private:
    FloatType inv_count_;
    FloatType inv_countm1_;
    af::versa<FloatType, af::c_grid<2> > sum_;
    af::versa<FloatType, af::c_grid<2> > sq_sum_;
  };

  /**
   * Calculate the masked mean and variance box filtered image. The
   * minimum number of counts specifies which pixels near a masked region
   * will be used.
   *
   * If the minimum counts == 1 then all non-masked pixels will be used.
   * If the minimum counts == size of the kernel then only though pixels
   * where all surrounding pixels are on will be used etc.
   *
   * A min_count value <= 0 will be interpreted as the size of the kernel.
   *
   * If min_count is set to 1 then a call to sample_variance will result
   * in an exception.
   */
  template <typename FloatType = double>
  class MeanAndVarianceFilterMasked {
  public:
    typedef FloatType value_type;

    /**
     * Initialise the algorithm
     * @param image The image to filter
     * @param mask The mask to use (0 = off, 1 == on)
     * @param size The size of the filter kernel (2 * size + 1)
     * @param min_count The minimum counts to use
     */
    MeanAndVarianceFilterMasked(const af::const_ref<FloatType, af::c_grid<2> > &image,
                                const af::const_ref<int, af::c_grid<2> > &mask,
                                int2 size,
                                int min_count)
        : min_count_(min_count), mask_(mask.accessor()) {
      const FloatType BIG = (1 << 24);  // About 1.6m counts

      // Check the input is valid
      DIALS_ASSERT(size.all_gt(0));
      DIALS_ASSERT(image.accessor().all_gt(0));
      DIALS_ASSERT(image.accessor().all_eq(mask.accessor()));

      // Copy the mask array
      for (std::size_t i = 0; i < mask.size(); ++i) {
        mask_[i] = mask[i] && image[i] < BIG;
      }

      // Ensure the min counts are valid
      if (min_count_ <= 0) {
        min_count_ = (2 * size[0] + 1) * (2 * size[1] + 1);
      } else {
        DIALS_ASSERT(min_count_ <= (2 * size[0] + 1) * (2 * size[1] + 1)
                     && min_count_ > 1);
      }

      // Calculate the summed area under the mask
      summed_mask_ = summed_area<int>(mask, size);

      // Ensure that all masked pixels are zero in the image and update the mask
      af::versa<FloatType, af::c_grid<2> > temp(image.accessor(),
                                                af::init_functor_null<FloatType>());
      af::versa<FloatType, af::c_grid<2> > image_sq(image.accessor(),
                                                    af::init_functor_null<FloatType>());
      for (std::size_t i = 0; i < image.size(); ++i) {
        temp[i] = image[i] * (mask[i] != 0);
        mask_[i] *= (summed_mask_[i] >= min_count_);
        image_sq[i] = temp[i] * temp[i];
      }

      // Calculate the summed image and the summed image**2
      summed_image_ = summed_area<FloatType>(temp.const_ref(), size);
      summed_image_sq_ = summed_area<FloatType>(image_sq.const_ref(), size);
    }

    /**
     * @returns The mean filtered image
     */
    af::versa<FloatType, af::c_grid<2> > mean() const {
      af::versa<FloatType, af::c_grid<2> > m(summed_image_.accessor(), 0);
      for (std::size_t i = 0; i < summed_image_.size(); ++i) {
        if (mask_[i]) {
          m[i] = summed_image_[i] / summed_mask_[i];
        }
      }
      return m;
    }

    /**
     * @returns The variance filtered image
     */
    af::versa<FloatType, af::c_grid<2> > variance() const {
      af::versa<FloatType, af::c_grid<2> > v(summed_image_.accessor(), 0);
      for (std::size_t i = 0; i < summed_image_.size(); ++i) {
        if (mask_[i]) {
          int c = summed_mask_[i];
          double s = summed_image_[i];
          double s2 = summed_image_sq_[i];
          v[i] = (s2 - (s * s / c)) / (c);
        }
      }
      return v;
    }

    /**
     * @returns The sample variance filtered image.
     */
    af::versa<FloatType, af::c_grid<2> > sample_variance() const {
      DIALS_ASSERT(min_count_ > 1);
      af::versa<FloatType, af::c_grid<2> > v(summed_image_.accessor(), 0);
      for (std::size_t i = 0; i < summed_image_.size(); ++i) {
        if (mask_[i]) {
          int c = summed_mask_[i];
          FloatType s = summed_image_[i];
          FloatType s2 = summed_image_sq_[i];
          v[i] = (s2 - (s * s / c)) / (c - 1);
        }
      }
      return v;
    }

    /**
     * @returns The counts per local area
     */
    af::versa<int, af::c_grid<2> > count() const {
      return summed_mask_;
    }

    /**
     * @returns The update mask
     */
    af::versa<int, af::c_grid<2> > mask() const {
      return mask_;
    }

  private:
    int min_count_;
    af::versa<int, af::c_grid<2> > mask_;
    af::versa<int, af::c_grid<2> > summed_mask_;
    af::versa<FloatType, af::c_grid<2> > summed_image_;
    af::versa<FloatType, af::c_grid<2> > summed_image_sq_;
  };

}}  // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_IMAGE_FILTER_MEAN_AND_VARIANCE_H */
