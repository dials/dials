/*
 * normal_discriminator.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_BACKGROUND_NORMAL_DISCRIMINATOR_H
#define DIALS_ALGORITHMS_BACKGROUND_NORMAL_DISCRIMINATOR_H

#include <algorithm>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/ref_reductions.h>
#include <boost/math/special_functions/erf.hpp>
#include <scitbx/math/mean_and_variance.h>
#include <dials/array_family/sort_index.h>
#include <dials/error.h>
#include "discriminator_strategy.h"

namespace dials { namespace algorithms {

  using boost::math::erf_inv;
  using scitbx::af::shared;
  using scitbx::af::const_ref;
  using scitbx::af::ref;
  using scitbx::af::min;
  using scitbx::af::max;
  using scitbx::af::mean;
  using scitbx::math::mean_and_variance;
  using scitbx::af::flex_int;
  using scitbx::af::flex_double;
  using dials::af::sort_index;

  /**
   * Get the expected number of standard deviations based on the number of
   * observations. Given by erf(x / sqrt(2)) = 1 - 1 / N
   * @param n_obs The number of observations
   * @returns The expected number of standard deviations
   */
  inline
  double normal_expected_n_sigma(int n_obs) {
    return sqrt(2.0) * erf_inv(1.0 - (1.0 / n_obs));
  }

  /**
   * Get the maximum number of standard deviations in the range of data
   * @param n_obs The number of observations
   * @returns The expected number of standard deviations
   */
  double minimum_n_sigma(const const_ref<double> &data) {

    // Calculate the mean and standard deviation of the data
    mean_and_variance <double> mean_and_variance(data);
    double mean = mean_and_variance.mean();
    double sdev = mean_and_variance.unweighted_sample_standard_deviation();

    // If sdev is zero then the extent of the data is 0 sigma
    if (sdev == 0) {
      return 0.0;
    }

    // Calculate the min/max of the data
    double mind = min(data);

    // Calculate t-statistic of min/max
    double min_n_sigma = (mind - mean) / sdev;

    // return the maximum number of sigma
    return min_n_sigma;
  }

  /**
   * Get the maximum number of standard deviations in the range of data
   * @param n_obs The number of observations
   * @returns The expected number of standard deviations
   */
  double maximum_n_sigma(const const_ref<double> &data) {

    // Calculate the mean and standard deviation of the data
    mean_and_variance <double> mean_and_variance(data);
    double mean = mean_and_variance.mean();
    double sdev = mean_and_variance.unweighted_sample_standard_deviation();

    // If sdev is zero then the extent of the data is 0 sigma
    if (sdev == 0) {
      return 0.0;
    }

    // Calculate the min/max of the data
    double maxd = max(data);

    // Calculate t-statistic of min/max
    double max_n_sigma = (maxd - mean) / sdev;

    // return the maximum number of sigma
    return max_n_sigma;
  }

  /**
   * Get the maximum number of standard deviations in the range of data
   * @param n_obs The number of observations
   * @returns The expected number of standard deviations
   */
  double absolute_maximum_n_sigma(const const_ref<double> &data) {

    // Calculate the mean and standard deviation of the data
    mean_and_variance <double> mean_and_variance(data);
    double mean = mean_and_variance.mean();
    double sdev = mean_and_variance.unweighted_sample_standard_deviation();

    // If sdev is zero then the extent of the data is 0 sigma
    if (sdev == 0) {
      return 0.0;
    }

    // Calculate the min/max of the data
    double mind = min(data);
    double maxd = max(data);

    // Calculate t-statistic of min/max
    double min_n_sigma = (mean - mind) / sdev;
    double max_n_sigma = (maxd - mean) / sdev;

    //std::cout << min_n_sigma << " " << max_n_sigma << std::endl;

    // return the maximum number of sigma
    return max_n_sigma > min_n_sigma ? max_n_sigma : min_n_sigma;
  }

  /**
   * Check if the data is normally distributed.
   *
   * Calculate the t-statistic of the min/max of the data and check if it is
   * between the given n_sigma.
   *
   * @param data The array of pixel values
   * @returns True/False
   */
  inline
  bool is_normally_distributed(const const_ref<double> &data, double n_sigma) {

    // Get the maximum n sigma
    double max_n_sigma = absolute_maximum_n_sigma(data);

    // return whether within required sigma
    return max_n_sigma < n_sigma;
  }

  /**
   * Check if the data is normally distributed.
   *
   * Calculate the t-statistic of the min/max of the data and check if it is
   * between the expected n_sigma
   *
   * @param data The array of pixel values
   * @returns True/False
   */
  inline
  bool is_normally_distributed(const const_ref<double> &data) {

    // Calculate expected sigma from number of points
    double n_sigma = normal_expected_n_sigma(data.size());

    // Check if data is normally distributed using sigma value
    return is_normally_distributed(data, n_sigma);
  }


  /**
   * A class that uses normal distribution statistics to discriminate
   * between background and peak pixels in a reflection shoebox.
   */
  class NormalDiscriminator : public DiscriminatorStrategy {
  public:

    /** Initialise the class. */
    NormalDiscriminator()
      : min_data_(10),
        n_sigma_(3.0) {}

    /** Initialise the class with parameters. */
    NormalDiscriminator(std::size_t min_data, double n_sigma)
      : min_data_(min_data),
        n_sigma_(n_sigma) {
      DIALS_ASSERT(min_data > 0);
      DIALS_ASSERT(n_sigma > 0.0);
    }

    /**
     * Discriminate between peak and background pixels.
     *
     * First get the indices of those pixels that belong to the reflection.
     * Sort the pixels in order of ascending intensity. Then check if the
     * intensities are normally distributed. If not then remove the pixel
     * with the highest intensity from the list and check again. Keep going
     * untill the list of pixels is normally distributed, or the maximum
     * number of iterations is reached. The remaining pixels are classes
     * as background, the rest are peak.
     *
     * @params shoebox The shoebox profile
     * @params mask The shoebox mask
     */
    void operator()(const flex_int &shoebox, flex_int &mask) const {

      // Ensure data is correctly sized.
      DIALS_ASSERT(shoebox.size() == mask.size());

      // Copy valid pixels and indices into list
      shared<int> indices;
      for (std::size_t i = 0; i < shoebox.size(); ++i) {
        if (mask[i]) {
          indices.push_back(i);
        }
      }

      // Check we have enough data
      DIALS_ASSERT(indices.size() >= min_data_);

      // Sort the pixels into ascending intensity order
      sort_index(indices.begin(), indices.end(), shoebox.begin());
      flex_double pixels(indices.size());
      for (std::size_t i = 0; i < indices.size(); ++i) {
        pixels[i] = (double)shoebox[indices[i]];
      }

      // Check if the data is normally distributed. If it is not, then remove
      // a value of high intensity and keep looping until it is. If the number
      // of iterations exceeds the maximum then exit the loop.
      std::size_t num_data = pixels.size();
      for (; num_data > min_data_; --num_data) {
        if (is_normally_distributed(const_ref<double>(
            pixels.begin(), num_data), n_sigma_)) {
          break;
        }
      }

      // Set all the rejected pixels as peak pixels and all the accepted
      // pixels as background pixels
      for (std::size_t i = 0; i < num_data; ++i) {
        mask[indices[i]] = (1 << 0);
      }
      for (std::size_t i = num_data; i < indices.size(); ++i) {
        mask[indices[i]] = (1 << 1);
      }
    }

    /**
     * Process just a shoebox and return a mask
     * @param shoebox The shoebox profile
     * @return The mask
     */
    flex_int operator()(const flex_int &shoebox) const {
      flex_int mask(shoebox.accessor(), 1);
      this->operator()(shoebox, mask);
      return mask;
    }

    /**
     * Process the reflection
     * @param reflection The reflection
     */
    virtual void operator()(Reflection &reflection) const {
      this->operator()(reflection.get_shoebox(), reflection.get_shoebox_mask());
    }

  private:

    std::size_t min_data_;
    double n_sigma_;
  };
}}

#endif /* DIALS_ALGORITHMS_BACKGROUND_NORMAL_DISCRIMINATOR_H */
