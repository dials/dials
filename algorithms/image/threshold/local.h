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
  inline
  af::versa< bool, af::c_grid<2> > niblack(
      const af::const_ref< double, af::c_grid<2> > &image,
      int2 size, double n_sigma) {

    // Check the input
    DIALS_ASSERT(n_sigma >= 0);

    // Calculate the mean and variance filtered images
    MeanAndVarianceFilter filter(image, size);
    af::versa< double, af::c_grid<2> > mean = filter.mean();
    af::versa< double, af::c_grid<2> > var = filter.sample_variance();

    // Assign the pixels to object and background
    af::versa< bool, af::c_grid<2> > result(image.accessor());
    #pragma omp parallel for
    for (std::size_t i = 0; i < var.size(); ++i) {
      result[i] = image[i] > mean[i] + n_sigma * sqrt(var[i]) ? 1 : 0;
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
  inline
  af::versa< bool, af::c_grid<2> > sauvola(
      const af::const_ref< double, af::c_grid<2> > &image,
      int2 size, double k, double r) {
    // Check the input
    DIALS_ASSERT(k >= 0 && r >= 0);

    // Calculate the mean and variance filtered images
    MeanAndVarianceFilter filter(image, size);
    af::versa< double, af::c_grid<2> > mean = filter.mean();
    af::versa< double, af::c_grid<2> > var = filter.sample_variance();

    // Assign the pixels to object and background
    af::versa< bool, af::c_grid<2> > result(image.accessor());
    #pragma omp parallel for
    for (std::size_t i = 0; i < var.size(); ++i) {
      result[i] = image[i] > mean[i] * (
        1.0 + k * (sqrt(var[i]) / r - 1)) ? 1 : 0;
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
  inline
  af::versa< bool, af::c_grid<2> > fano(
      const af::const_ref< double, af::c_grid<2> > &image,
      int2 size, double n_sigma) {

    // Check the input
    DIALS_ASSERT(n_sigma >= 0);

    // Calculate the fano filtered image
    FanoFilter filter(image, size);
    af::versa< double, af::c_grid<2> > fano_image = filter.fano();
    af::versa< int, af::c_grid<2> >    fano_mask  = filter.mask();

    // Calculate the bound
    std::size_t n = (2 * size[0] + 1) * (2 * size[1] + 1);
    double bound = 1.0 + n_sigma * sqrt(2.0 / (n - 1));

    // Assign pixels to object or background
    af::versa< bool, af::c_grid<2> > result(image.accessor());
    #pragma omp parallel for
    for (std::size_t i = 0; i < image.size(); ++i) {
      result[i] = (fano_mask[i] && fano_image[i]) > bound ? 1 : 0;
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
  inline
  af::versa< bool, af::c_grid<2> > fano_masked(
      const af::const_ref< double, af::c_grid<2> > &image,
      const af::const_ref< bool, af::c_grid<2> > &mask,
      int2 size, int min_count, double n_sigma) {

    // Check the input
    DIALS_ASSERT(n_sigma >= 0);

    // Copy the mask into a temp variable
    af::versa< int, af::c_grid<2> > temp(mask.accessor());
    #pragma omp parallel for
    for (std::size_t i = 0; i < temp.size(); ++i) {
      temp[i] = mask[i] ? 1 : 0;
    }

    // Calculate the masked fano filtered image
    FanoFilterMasked filter(image, temp.const_ref(), size, min_count);
    af::versa< double, af::c_grid<2> > fano_image = filter.fano();
    af::versa< int, af::c_grid<2> >    count      = filter.count();
    temp                                          = filter.mask();

    // Assign pixels to object or background
    af::versa< bool, af::c_grid<2> > result(image.accessor(), false);
    #pragma omp parallel for
    for (std::size_t i = 0; i < image.size(); ++i) {
      if (temp[i]) {
        double bound = 1.0 + n_sigma * sqrt(2.0 / (count[i] - 1));
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
  inline
  af::versa< bool, af::c_grid<2> > gain(
      const af::const_ref< double, af::c_grid<2> > &image,
      const af::const_ref< bool, af::c_grid<2> > &mask,
      af::ref< double, af::c_grid<2> > gain, int2 size,
      int min_count, double n_sigma) {

    // Check the input
    DIALS_ASSERT(n_sigma >= 0);

    // Copy the mask into a temp variable
    af::versa< int, af::c_grid<2> > temp(mask.accessor());
    #pragma omp parallel for
    for (std::size_t i = 0; i < temp.size(); ++i) {
      temp[i] = mask[i] ? 1 : 0;
    }

    // Calculate the masked fano filtered image
    FanoFilterMasked filter(image, temp.const_ref(), size, min_count);
    af::versa< double, af::c_grid<2> > fano_image = filter.fano();
    af::versa< int, af::c_grid<2> >    count      = filter.count();
    temp                                          = filter.mask();

    // Assign pixels to object or background
    af::versa< bool, af::c_grid<2> > result(image.accessor(), false);
    #pragma omp parallel for
    for (std::size_t i = 0; i < image.size(); ++i) {
      if (temp[i]) {
        double bound = gain[i] + n_sigma * gain[i] * sqrt(2.0 / (count[i] - 1));
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
   */
  inline
  af::versa< bool, af::c_grid<2> > kabsch(
      const af::const_ref< double, af::c_grid<2> > &image,
      const af::const_ref< bool, af::c_grid<2> > &mask,
      int2 size, double nsig_b, double nsig_s) {

    // Check the input
    DIALS_ASSERT(nsig_b >= 0 && nsig_s >= 0);

    // Copy the mask into a temp variable
    af::versa< int, af::c_grid<2> > temp(mask.accessor());
    #pragma omp parallel for
    for (std::size_t i = 0; i < temp.size(); ++i) {
      temp[i] = mask[i] ? 1 : 0;
    }

    // Calculate the masked fano filtered image
    FanoFilterMasked filter(image, temp.const_ref(), size, 0);
    af::versa< double, af::c_grid<2> > fano_image = filter.fano();
    af::versa< double, af::c_grid<2> > mean       = filter.mean();
    af::versa< int, af::c_grid<2> > count         = filter.count();
    temp                                          = filter.mask();

    // Assign pixels to object or background
    af::versa< bool, af::c_grid<2> > result(image.accessor(), 0);
    #pragma omp parallel for
    for (std::size_t i = 0; i < image.size(); ++i) {
      if (temp[i]) {
        double bnd_b = 1.0 + nsig_b * sqrt(2.0 / (count[i] - 1));
        double bnd_s = mean[i] + nsig_s * sqrt(mean[i]);
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
   */
  inline
  af::versa< bool, af::c_grid<2> > kabsch_w_gain(
      const af::const_ref< double, af::c_grid<2> > &image,
      const af::const_ref< bool, af::c_grid<2> > &mask,
      af::ref< double, af::c_grid<2> > gain, int2 size,
      double nsig_b, double nsig_s) {

    // Check the input
    DIALS_ASSERT(nsig_b >= 0 && nsig_s >= 0);

    // Copy the mask into a temp variable
    af::versa< int, af::c_grid<2> > temp(mask.accessor());
    #pragma omp parallel for
    for (std::size_t i = 0; i < temp.size(); ++i) {
      temp[i] = mask[i] ? 1 : 0;
    }

    // Calculate the masked fano filtered image
    FanoFilterMasked filter(image, temp.const_ref(), size, 0);
    af::versa< double, af::c_grid<2> > fano_image = filter.fano();
    af::versa< double, af::c_grid<2> > mean       = filter.mean();
    af::versa< int, af::c_grid<2> > count         = filter.count();
    temp                                          = filter.mask();

    // Assign pixels to object or background
    af::versa< bool, af::c_grid<2> > result(image.accessor(), 0);
    #pragma omp parallel for
    for (std::size_t i = 0; i < image.size(); ++i) {
      if (temp[i]) {
        double bnd_b = gain[i] + nsig_b * gain[i] * sqrt(2.0 / (count[i] - 1));
        double bnd_s = mean[i] + nsig_s * sqrt(gain[i] * mean[i]);
        result[i]  = (fano_image[i] > bnd_b && image[i] > bnd_s) ? 1 : 0;
      }
    }

    // Return thresholded image
    return result;
  }

}} // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_IMAGE_THRESHOLD_UNIMODAL_H */
