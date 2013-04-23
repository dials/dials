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

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/ref_reductions.h>
#include <dials/error.h>
#include "summed_area.h"

namespace dials { namespace algorithms {

  using scitbx::af::int2;
  using scitbx::af::flex;
  using scitbx::af::flex_int;
  using scitbx::af::flex_bool;
  using scitbx::af::flex_double;

  /**
   * Calculate the mean box filtered image.
   * @param image The image to filter
   * @param size The size of the filter kernel (2 * size + 1)
   * @returns The filtered image
   */
  inline
  flex_double mean_filter(const flex_double &image, int2 size) {
    flex_double mean = summed_area<double>(image, size);
    double inv_count = 1.0 / ((double)(size[0] * size[1]));
    
    #pragma omp parallel for
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
   * @returns The filtered image
   */
  inline
  flex_double mean_filter_masked(const flex_double &image,
      flex_int &mask, int2 size, int min_count) {

    // Check the input is valid
    DIALS_ASSERT(size[0] > 0 && size[1] > 0);
    DIALS_ASSERT(image.accessor().all().size() == 2);
    DIALS_ASSERT(image.accessor().all().all_gt(0));
    DIALS_ASSERT(image.accessor().all().all_eq(mask.accessor().all()));

    // Ensure the min counts are valid
    if (min_count <= 0) {
      min_count = (2 * size[0] + 1) * (2 * size[1] + 1);
    } else {
      DIALS_ASSERT(min_count <= (2 * size[0] + 1) * (2 * size[1] + 1));
    }

    // Calculate the summed area under the mask
    flex_int summed_mask = summed_area<int>(mask, size);

    // Ensure that all masked pixels are zero in the image and update the mask
    flex_double temp(image);
    #pragma omp parallel for    
    for (std::size_t i = 0; i < image.size(); ++i) {
      temp[i] *= (mask[i] != 0);
      mask[i] *= (summed_mask[i] >= min_count);
    }

    // Calculate the summed area under the image
    flex_double summed_image = summed_area<double>(temp, size);

    // Calculate the mean filtered image
    #pragma omp parallel for    
    for (std::size_t i = 0; i < image.size(); ++i) {
      if (mask[i]) {
        summed_image[i] /= summed_mask[i];
      } else {
        summed_image[i] = 0;
      }
    }

    // Return the filtered image
    return summed_image;
  }


  /**
   * A class to calculate the mean and variance filtered images.
   */
  class MeanAndVarianceFilter {
  public:

    /**
     * Initialise the algorithm.
     * @params image The image to filter
     * @param size The size of the filter kernel (2 * size + 1)
     */
    MeanAndVarianceFilter(const flex_double &image, int2 size) {

      // Check the input is valid
      DIALS_ASSERT(size[0] > 0 && size[1] > 0);
      DIALS_ASSERT(image.accessor().all().size() == 2);
      DIALS_ASSERT(image.accessor().all().all_gt(0));

      // Calculate the squared image
      flex_double image_sq(image.accessor());
      #pragma omp parallel for      
      for (std::size_t i = 0; i < image.size(); ++i) {
        image_sq[i] = image[i] * image[i];
      }

      // Inverse of counts to avoid excessive division
      int count = (size[0] * size[1]);
      inv_count_ = 1.0 / count;
      inv_countm1_ = 1.0 / (count - 1);

      // Calculate the summed area under the image and image**2
      sum_ = summed_area<double>(image, size);
      sq_sum_ = summed_area<double>(image_sq, size);
    }

    /**
     * Calculate the mean filtered image
     * @returns The mean filtered image.
     */
    flex_double mean() const {
      flex_double m(sum_.accessor());
      #pragma omp parallel for        
      for (std::size_t i = 0; i < sum_.size(); ++i) {
        m[i] = sum_[i] * inv_count_;
      }
      return m;
    }

    /**
     * Calculate the variance filtered image
     * @returns The variance filtered image
     */
    flex_double variance() const {
      flex_double v(sum_.accessor());
      #pragma omp parallel for        
      for (std::size_t i = 0; i < sum_.size(); ++i) {
        v[i] = (sq_sum_[i] - (sum_[i] * sum_[i] * inv_count_)) * (inv_count_);
      }
      return v;
    }

    /**
     * Calculate the sample variance filtered image
     * @returns The sample variance filtered image
     */
    flex_double sample_variance() const {
      flex_double v(sum_.accessor());
      #pragma omp parallel for        
      for (std::size_t i = 0; i < sum_.size(); ++i) {
        v[i] = (sq_sum_[i] - (sum_[i] * sum_[i] * inv_count_)) * (inv_countm1_);
      }
      return v;
    }

  private:
    double inv_count_;
    double inv_countm1_;
    flex_double sum_;
    flex_double sq_sum_;
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
  class MeanAndVarianceFilterMasked {
  public:

    /**
     * Initialise the algorithm
     * @param image The image to filter
     * @param mask The mask to use (0 = off, 1 == on)
     * @param size The size of the filter kernel (2 * size + 1)
     * @param min_count The minimum counts to use
     */
    MeanAndVarianceFilterMasked(const flex_double &image,
        const flex_int &mask, int2 size, int min_count)
        : min_count_(min_count), mask_(mask) {

      // Check the input is valid
      DIALS_ASSERT(size[0] > 0 && size[1] > 0);
      DIALS_ASSERT(image.accessor().all().size() == 2);
      DIALS_ASSERT(image.accessor().all().all_gt(0));
      DIALS_ASSERT(image.accessor().all().all_eq(mask.accessor().all()));

      // Ensure the min counts are valid
      if (min_count_ <= 0) {
        min_count_ = (2 * size[0] + 1) * (2 * size[1] + 1);
      } else {
        DIALS_ASSERT(min_count_ <= (2 * size[0] + 1) * (2 * size[1] + 1));
      }

      // Calculate the summed area under the mask
      summed_mask_ = summed_area<int>(mask, size);

      // Ensure that all masked pixels are zero in the image and update the mask
      flex_double temp(image);
      flex_double image_sq(image.accessor());
      #pragma omp parallel for        
      for (std::size_t i = 0; i < image.size(); ++i) {
        temp[i] *= (mask[i] != 0);
        mask_[i] *= (summed_mask_[i] >= min_count_);
        image_sq[i] = temp[i] * temp[i];
      }

      // Calculate the summed image and the summed image**2
      summed_image_ = summed_area<double>(temp, size);
      summed_image_sq_ = summed_area<double>(image_sq, size);
    }

    /**
     * @returns The mean filtered image
     */
    flex_double mean() const {
      flex_double m(summed_image_.accessor(), 0);
      #pragma omp parallel for        
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
    flex_double variance() const {
      flex_double v(summed_image_.accessor(), 0);
      #pragma omp parallel for        
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
    flex_double sample_variance() const {
      DIALS_ASSERT(min_count_ > 1);
      flex_double v(summed_image_.accessor(), 0);
      #pragma omp parallel for        
      for (std::size_t i = 0; i < summed_image_.size(); ++i) {
        if (mask_[i]) {
          int c = summed_mask_[i];
          double s = summed_image_[i];
          double s2 = summed_image_sq_[i];
          v[i] = (s2 - (s * s / c)) / (c - 1);
        }
      }
      return v;
    }

    /**
     * @returns The counts per local area
     */
    flex_int count() const {
      return summed_mask_;
    }

    /**
     * @returns The update mask
     */
    flex_int mask() const {
      return mask_;
    }

  private:
    int min_count_;
    flex_int mask_;
    flex_int summed_mask_;
    flex_double summed_image_;
    flex_double summed_image_sq_;
  };

}} // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_IMAGE_FILTER_MEAN_AND_VARIANCE_H */
