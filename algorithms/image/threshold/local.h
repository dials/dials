/*
 * local.h
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

#include <omptbx/omp_or_stubs.h>

#include <cmath>
#include <iostream>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/ref_reductions.h>
#include <dials/error.h>
#include <dials/algorithms/image/filter/mean_and_variance.h>
#include <dials/algorithms/image/filter/fano_filter.h>

namespace dials { namespace algorithms {

  /**
   * Threshold the image using the niblack method.
   *
   * pixel > mean + n_sigma * sdev ? object : background
   *
   * @param image The image to threshold
   * @param size The size of the local area
   * @param n_sigma The number of standard deviations
   * @returns The thresholded image
   */
  template <typename FloatType>
  af::versa< bool, af::c_grid<2> > niblack(
      const af::const_ref< FloatType, af::c_grid<2> > &image,
      int2 size, double n_sigma) {

    // Check the input
    DIALS_ASSERT(n_sigma >= 0);

    // Calculate the mean and variance filtered images
    MeanAndVarianceFilter<FloatType> filter(image, size);
    af::versa< FloatType, af::c_grid<2> > mean = filter.mean();
    af::versa< FloatType, af::c_grid<2> > var = filter.sample_variance();

    // Assign the pixels to object and background
    af::versa< bool, af::c_grid<2> > result(image.accessor(),
      af::init_functor_null<bool>());
    #pragma omp parallel for
    for (std::size_t i = 0; i < var.size(); ++i) {
      result[i] = image[i] > mean[i] + n_sigma * std::sqrt(var[i]) ? 1 : 0;
    }

    // Return the thresholded image
    return result;
  }

  /**
   * Threshold the image using the sauvola method.
   *
   * pixel > mean * (1 + k * (sdev / (r - 1)))
   *
   * @param image The image to threshold
   * @param size The size of the local area
   * @param k A parameter
   * @param r A parameter
   * @returns The thresholded image
   */
  template <typename FloatType>
  af::versa< bool, af::c_grid<2> > sauvola(
      const af::const_ref< FloatType, af::c_grid<2> > &image,
      int2 size, double k, double r) {
    // Check the input
    DIALS_ASSERT(k >= 0 && r > 1);

    // Calculate the mean and variance filtered images
    MeanAndVarianceFilter<FloatType> filter(image, size);
    af::versa< FloatType, af::c_grid<2> > mean = filter.mean();
    af::versa< FloatType, af::c_grid<2> > var = filter.sample_variance();

    // Assign the pixels to object and background
    af::versa< bool, af::c_grid<2> > result(image.accessor(),
      af::init_functor_null<bool>());
    #pragma omp parallel for
    for (std::size_t i = 0; i < var.size(); ++i) {
      result[i] = image[i] > mean[i] * (
        1.0 + k * (std::sqrt(var[i]) / r - 1)) ? 1 : 0;
    }

    // Return the thresholded image
    return result;
  }


  /**
   * Threshold the image using a fano filter. Essentially a test for objects
   * within a poisson distribution.
   *
   * var/mean > 1.0 + n_sigma * sqrt(2 / (n - 1)) ? object : background
   *
   * @param image The image to threshold
   * @param size The size of the local window
   * @param n_sigma The number of standard deviations.
   */
  template <typename FloatType>
  af::versa< bool, af::c_grid<2> > fano(
      const af::const_ref< FloatType, af::c_grid<2> > &image,
      int2 size, double n_sigma) {

    // Check the input
    DIALS_ASSERT(n_sigma >= 0);

    // Calculate the fano filtered image
    FanoFilter<FloatType> filter(image, size);
    af::versa< FloatType, af::c_grid<2> > fano_image = filter.fano();

    // Calculate the bound
    std::size_t n = (2 * size[0] + 1) * (2 * size[1] + 1);
    DIALS_ASSERT(n > 1);
    FloatType bound = 1.0 + n_sigma * std::sqrt(2.0 / (n - 1));

    // Assign pixels to object or background
    af::versa< bool, af::c_grid<2> > result(image.accessor(),
      af::init_functor_null<bool>());
    #pragma omp parallel for
    for (std::size_t i = 0; i < image.size(); ++i) {
      result[i] = (fano_image[i] > bound) ? 1 : 0;
    }

    // Return thresholded image
    return result;
  }

  /**
   * Threshold the image using a fano filter. Essentially a test for objects
   * within a poisson distribution.
   *
   * var/mean > 1.0 + n_sigma * sqrt(2 / (n - 1)) ? object : background
   *
   * @param image The image to threshold
   * @param mask The mask to use
   * @param size The size of the local window
   * @param min_count The minimum counts for a point to be valid
   * @param n_sigma The number of standard deviations.
   */
  template <typename FloatType>
  af::versa< bool, af::c_grid<2> > fano_masked(
      const af::const_ref< FloatType, af::c_grid<2> > &image,
      const af::const_ref< bool, af::c_grid<2> > &mask,
      int2 size, int min_count, double n_sigma) {

    // Check the input
    DIALS_ASSERT(n_sigma >= 0);
    DIALS_ASSERT(min_count > 1);

    // Copy the mask into a temp variable
    af::versa< int, af::c_grid<2> > temp(mask.accessor());
    #pragma omp parallel for
    for (std::size_t i = 0; i < temp.size(); ++i) {
      temp[i] = mask[i] ? 1 : 0;
    }

    // Calculate the masked fano filtered image
    FanoFilterMasked<FloatType> filter(image, temp.const_ref(), size, min_count);
    af::versa< FloatType, af::c_grid<2> > fano_image = filter.fano();
    af::versa< int, af::c_grid<2> > count = filter.count();
    temp = filter.mask();

    // Assign pixels to object or background
    af::versa< bool, af::c_grid<2> > result(image.accessor(), false);
    #pragma omp parallel for
    for (std::size_t i = 0; i < image.size(); ++i) {
      if (temp[i]) {
        FloatType bound = 1.0 + n_sigma * std::sqrt(2.0 / (count[i] - 1));
        result[i] = (fano_image[i] > bound) ? 1 : 0;
      }
    }

    // Return thresholded image
    return result;
  }

  /**
   * Threshold the image using a gain filter. Same as the fano filter but
   * using a gain map for the calculation
   *
   * var/mean > g + n_sigma * g * sqrt(2 / (n - 1)) ? object : background
   *
   * @param image The image to threshold
   * @param mask The mask to use
   * @param gain The gain map
   * @param size The size of the local window
   * @param min_count The minimum counts for a point to be valid
   * @param n_sigma The number of standard deviations.
   */
  template <typename FloatType>
  af::versa< bool, af::c_grid<2> > gain(
      const af::const_ref< FloatType, af::c_grid<2> > &image,
      const af::const_ref< bool, af::c_grid<2> > &mask,
      const af::const_ref< FloatType, af::c_grid<2> > &gain,
      int2 size, int min_count, double n_sigma) {

    // Check the input
    DIALS_ASSERT(n_sigma >= 0);
    DIALS_ASSERT(min_count > 1);

    // Copy the mask into a temp variable
    af::versa< int, af::c_grid<2> > temp(mask.accessor());
    #pragma omp parallel for
    for (std::size_t i = 0; i < temp.size(); ++i) {
      temp[i] = mask[i] ? 1 : 0;
    }

    // Calculate the masked fano filtered image
    FanoFilterMasked<FloatType> filter(image, temp.const_ref(), size, min_count);
    af::versa< FloatType, af::c_grid<2> > fano_image = filter.fano();
    af::versa< int, af::c_grid<2> > count = filter.count();
    temp = filter.mask();

    // Assign pixels to object or background
    af::versa< bool, af::c_grid<2> > result(image.accessor(), false);
    #pragma omp parallel for
    for (std::size_t i = 0; i < image.size(); ++i) {
      if (temp[i]) {
        FloatType bound = gain[i] + n_sigma * gain[i] *
          std::sqrt(2.0 / (count[i] - 1));
        result[i] = (fano_image[i] > bound) ? 1 : 0;
      }
    }

    // Return thresholded image
    return result;
  }

  /**
   * Threshold the image as in xds. Same as the fano filter but
   * using a gain map for the calculation
   *
   * var/mean > g + n_sigma * g * sqrt(2 / (n - 1)) &&
   * pixel > mean + sqrt(mean) ? object : background
   *
   * @param image The image to threshold
   * @param mask The mask to use
   * @param size The size of the local window
   * @param nsig_b The background threshold.
   * @param nsig_s The strong pixel threshold
   * @param min_count The minimum number of pixels in the local area
   */
  template <typename FloatType>
  af::versa< bool, af::c_grid<2> > kabsch(
      const af::const_ref< FloatType, af::c_grid<2> > &image,
      const af::const_ref< bool, af::c_grid<2> > &mask,
      int2 size, double nsig_b, double nsig_s, int min_count) {

    // Check the input
    DIALS_ASSERT(nsig_b >= 0 && nsig_s >= 0);

    // Copy the mask into a temp variable
    af::versa< int, af::c_grid<2> > temp(mask.accessor());
    #pragma omp parallel for
    for (std::size_t i = 0; i < temp.size(); ++i) {
      temp[i] = mask[i] ? 1 : 0;
    }

    // Calculate the masked fano filtered image
    FanoFilterMasked<FloatType> filter(image, temp.const_ref(), size, min_count);
    af::versa< FloatType, af::c_grid<2> > fano_image = filter.fano();
    af::versa< FloatType, af::c_grid<2> > mean = filter.mean();
    af::versa< int, af::c_grid<2> > count = filter.count();
    temp = filter.mask();

    // Assign pixels to object or background
    af::versa< bool, af::c_grid<2> > result(image.accessor(), false);
    #pragma omp parallel for
    for (std::size_t i = 0; i < image.size(); ++i) {
      if (temp[i]) {
        FloatType bnd_b = 1.0 + nsig_b * std::sqrt(2.0 / (count[i] - 1));
        FloatType bnd_s = mean[i] + nsig_s * std::sqrt(mean[i]);
        result[i]  = (fano_image[i] > bnd_b && image[i] > bnd_s) ? 1 : 0;
      }
    }

    // Return thresholded image
    return result;
  }

  /**
   * Threshold the image as in xds. Same as the fano filter but
   * using a gain map for the calculation
   *
   * var/mean > g + n_sigma * g * sqrt(2 / (n - 1)) &&
   * pixel > mean + sqrt(gain * mean) ? object : background
   *
   * @param image The image to threshold
   * @param mask The mask to use
   * @param gain The gain map
   * @param size The size of the local window
   * @param nsig_b The background threshold.
   * @param nsig_s The strong pixel threshold
   * @param min_count The minimum number of pixels in the local area
   */
  template <typename FloatType>
  af::versa< bool, af::c_grid<2> > kabsch_w_gain(
      const af::const_ref< FloatType, af::c_grid<2> > &image,
      const af::const_ref< bool, af::c_grid<2> > &mask,
      const af::const_ref< FloatType, af::c_grid<2> > &gain,
      int2 size, double nsig_b, double nsig_s, int min_count) {

    // Check the input
    DIALS_ASSERT(nsig_b >= 0 && nsig_s >= 0);

    // Copy the mask into a temp variable
    af::versa< int, af::c_grid<2> > temp(mask.accessor());
    #pragma omp parallel for
    for (std::size_t i = 0; i < temp.size(); ++i) {
      temp[i] = mask[i] ? 1 : 0;
    }

    // Calculate the masked fano filtered image
    FanoFilterMasked<FloatType> filter(image, temp.const_ref(), size, min_count);
    af::versa< FloatType, af::c_grid<2> > fano_image = filter.fano();
    af::versa< FloatType, af::c_grid<2> > mean = filter.mean();
    af::versa< int, af::c_grid<2> > count = filter.count();
    temp = filter.mask();

    // Assign pixels to object or background
    af::versa< bool, af::c_grid<2> > result(image.accessor(), false);
    #pragma omp parallel for
    for (std::size_t i = 0; i < image.size(); ++i) {
      if (temp[i]) {
        FloatType bnd_b = gain[i] + nsig_b * gain[i] *
          std::sqrt(2.0 / (count[i] - 1));
        FloatType bnd_s = mean[i] + nsig_s * std::sqrt(gain[i] * mean[i]);
        result[i]  = (fano_image[i] > bnd_b && image[i] > bnd_s) ? 1 : 0;
      }
    }

    // Return thresholded image
    return result;
  }

  /**
   * A class to help debug spot finding by exposing the results of various bits
   * of processing.
   */
  class KabschDebug {
  public:

    /**
     * Do the processing.
     * @param image The image array
     * @param mask The mask array
     * @param size The size of the local window
     * @param nsig_b The background threshold.
     * @param nsig_s The strong pixel threshold
     * @param min_count The minimum number of pixels in the local area
     */
    KabschDebug(const af::const_ref<double, af::c_grid<2> > &image,
                const af::const_ref<bool, af::c_grid<2> > &mask,
                int2 size,
                double nsig_b,
                double nsig_s,
                int min_count) {

      // Check the input
      DIALS_ASSERT(nsig_b >= 0 && nsig_s >= 0);

      // Copy the mask into a temp variable
      af::versa< int, af::c_grid<2> > temp(mask.accessor());
      for (std::size_t i = 0; i < temp.size(); ++i) {
        temp[i] = mask[i] ? 1 : 0;
      }

      // Calculate the masked fano filtered image
      FanoFilterMasked<double> filter(image, temp.const_ref(), size, min_count);
      mean_ = filter.mean();
      variance_ = filter.sample_variance();
      cv_ = filter.fano();
      af::versa< int, af::c_grid<2> > count = filter.count();
      temp = filter.mask();

      // Assign pixels to object or background
      cv_mask_ = af::versa< bool, af::c_grid<2> >(image.accessor(), false);
      value_mask_ = af::versa< bool, af::c_grid<2> >(image.accessor(), false);
      final_mask_ = af::versa< bool, af::c_grid<2> >(image.accessor(), false);
      for (std::size_t i = 0; i < image.size(); ++i) {
        if (temp[i]) {
          double bnd_b = 1.0 + nsig_b * std::sqrt(2.0 / (count[i] - 1));
          double bnd_s = mean_[i] + nsig_s * std::sqrt(mean_[i]);
          cv_mask_[i] = cv_[i] > bnd_b;
          value_mask_[i] = image[i] > bnd_s;
          final_mask_[i] = cv_mask_[i] && value_mask_[i];
        }
      }
    }

    /** @returns The mean map */
    af::versa<double, af::c_grid<2> > mean() const {
      return mean_;
    }

    /** @returns The variance map. */
    af::versa<double, af::c_grid<2> > variance() const {
      return variance_;
    }

    /** @returns The coefficient of variation map */
    af::versa<double, af::c_grid<2> > coefficient_of_variation() const {
      return cv_;
    }

    /** @returns The thresholded coefficient of variation mask */
    af::versa<bool, af::c_grid<2> > cv_mask() const {
      return cv_mask_;
    }

    /** @returns The thresholded value mask */
    af::versa<bool, af::c_grid<2> > value_mask() const {
      return value_mask_;
    }

    /** @returns The final mask of strong pixels. */
    af::versa<bool, af::c_grid<2> > final_mask() const {
      return final_mask_;
    }

  private:

    af::versa< double, af::c_grid<2> > mean_;
    af::versa< double, af::c_grid<2> > variance_;
    af::versa< double, af::c_grid<2> > cv_;
    af::versa< bool, af::c_grid<2> > cv_mask_;
    af::versa< bool, af::c_grid<2> > value_mask_;
    af::versa< bool, af::c_grid<2> > final_mask_;
  };

}} // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_IMAGE_THRESHOLD_UNIMODAL_H */
