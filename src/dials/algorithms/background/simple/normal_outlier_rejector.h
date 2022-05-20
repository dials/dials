/*
 * normal_outlier_rejector.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_BACKGROUND_NORMAL_OUTLIER_REJECTOR_H
#define DIALS_ALGORITHMS_BACKGROUND_NORMAL_OUTLIER_REJECTOR_H

#include <algorithm>
#include <scitbx/array_family/ref_reductions.h>
#include <boost/math/special_functions/erf.hpp>
#include <scitbx/math/mean_and_variance.h>
#include <dials/array_family/sort_index.h>
#include <dials/algorithms/shoebox/mask_code.h>
#include <dials/error.h>

namespace dials { namespace algorithms { namespace background {

  using boost::math::erf_inv;
  using dials::af::sort_index;
  using scitbx::af::max;
  using scitbx::af::mean;
  using scitbx::af::min;
  using scitbx::math::mean_and_variance;

  /**
   * Get the expected number of standard deviations based on the number of
   * observations. Given by erf(n_sdev / sqrt(2)) = 1 - 1 / n_obs.
   * This function returns the value of sqrt(2) * erf-1(1 - 1 / n_obs)
   * @param n_obs The number of observations
   * @returns The expected number of standard deviations
   */
  inline double normal_expected_n_sigma(int n_obs) {
    return std::sqrt(2.0) * erf_inv(1.0 - (1.0 / n_obs));
  }

  /**
   * Get the maximum number of standard deviations in the range of data
   * @param data The observations
   * @returns The expected number of standard deviations
   */
  template <typename FloatType>
  FloatType minimum_n_sigma(const af::const_ref<FloatType> &data) {
    // Calculate the mean and standard deviation of the data
    mean_and_variance<FloatType> mean_and_variance(data);
    FloatType mean = mean_and_variance.mean();
    FloatType sdev = mean_and_variance.unweighted_sample_standard_deviation();

    // If sdev is zero then the extent of the data is 0 sigma
    if (sdev == 0) {
      return 0.0;
    }

    // Calculate t-statistic of the min of the data
    return (min(data) - mean) / sdev;
  }

  /**
   * Get the maximum number of standard deviations in the range of data
   * @param data The observations
   * @returns The maximum number of standard deviations
   */
  template <typename FloatType>
  FloatType maximum_n_sigma(const af::const_ref<FloatType> &data) {
    // Calculate the mean and standard deviation of the data
    mean_and_variance<FloatType> mean_and_variance(data);
    FloatType mean = mean_and_variance.mean();
    FloatType sdev = mean_and_variance.unweighted_sample_standard_deviation();

    // If sdev is zero then the extent of the data is 0 sigma
    if (sdev == 0) {
      return 0.0;
    }

    // Calculate t-statistic of the max of the data
    return (max(data) - mean) / sdev;
  }

  /**
   * Get the maximum number of standard deviations in the range of data
   * @param data The observations
   * @returns The maximum number of standard deviations
   */
  template <typename FloatType>
  FloatType absolute_maximum_n_sigma(const af::const_ref<FloatType> &data) {
    // Calculate the mean and standard deviation of the data
    mean_and_variance<FloatType> mean_and_variance(data);
    FloatType mean = mean_and_variance.mean();
    FloatType sdev = mean_and_variance.unweighted_sample_standard_deviation();

    // If sdev is zero then the extent of the data is 0 sigma
    if (sdev == 0) {
      return 0.0;
    }

    // Calculate t-statistic of min/max
    FloatType min_n_sigma = (mean - min(data)) / sdev;
    FloatType max_n_sigma = (max(data) - mean) / sdev;
    return max_n_sigma > min_n_sigma ? max_n_sigma : min_n_sigma;
  }

  /**
   * Check if the data is normally distributed.
   *
   * Calculate the t-statistic of the min/max of the data and check if it is
   * between the given n_sigma.
   *
   * @param data The array of pixel values
   * @param n_sigma The number of standard deviations you expect
   * @returns True/False
   */
  template <typename FloatType>
  bool is_normally_distributed(const af::const_ref<FloatType> &data, double n_sigma) {
    return maximum_n_sigma(data) < n_sigma;
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
  template <typename FloatType>
  bool is_normally_distributed(const af::const_ref<FloatType> &data) {
    return is_normally_distributed<FloatType>(data,
                                              normal_expected_n_sigma(data.size()));
  }

  /**
   * A class that uses normal distribution statistics to discriminate
   * between background and peak pixels in a reflection shoebox.
   */
  class NormalOutlierRejector : public OutlierRejector {
  public:
    /**
     * @param min_data The minimum number of data points to use.
     * @param n_sigma The number of standard deviations to check for
     */
    NormalOutlierRejector(std::size_t min_data) : min_data_(min_data) {
      DIALS_ASSERT(min_data > 0);
    }

    /**
     * Discriminate between peak and background pixels.
     *
     * First get the indices of those pixels that belong to the reflection.
     * Sort the pixels in order of ascending intensity. Then check if the
     * intensities are normally distributed. If not then remove the pixel
     * with the highest intensity from the list and check again. Keep going
     * until the list of pixels is normally distributed, or the maximum
     * number of iterations is reached. The remaining pixels are classed
     * as background, the rest are peak.
     *
     * The mask is used in both input and output. On input the mask is checked
     * for valid pixels. The discriminated peak/background pixels are then
     * written into the mask. The pixel values to include in the discrimination
     * procedure must be set to Valid | Background.
     *
     * @params shoebox The shoebox profile
     * @params mask The shoebox mask
     */
    virtual void mark(const af::const_ref<double, af::c_grid<3> > &shoebox,
                      af::ref<int, af::c_grid<3> > mask) const {
      // Ensure data is correctly sized.
      DIALS_ASSERT(shoebox.size() == mask.size());

      // Copy valid pixels and indices into list
      int mask_code = shoebox::Valid | shoebox::Background;
      af::shared<int> indices;
      for (std::size_t i = 0; i < shoebox.size(); ++i) {
        if ((mask[i] & mask_code) == mask_code
            && (mask[i] & shoebox::Overlapped) == 0) {
          indices.push_back(i);
        }
      }

      // Check we have enough data
      DIALS_ASSERT(indices.size() >= min_data_);

      // Sort the pixels into ascending intensity order
      sort_index(indices.begin(), indices.end(), shoebox.begin());
      af::shared<double> pixels(indices.size(), af::init_functor_null<double>());
      for (std::size_t i = 0; i < indices.size(); ++i) {
        pixels[i] = (double)shoebox[indices[i]];
      }

      // Check if the data is normally distributed. If it is not, then remove
      // a value of high intensity and keep looping until it is. If the number
      // of iterations exceeds the maximum then exit the loop.
      std::size_t num_data = pixels.size();
      for (; num_data > min_data_; --num_data) {
        if (is_normally_distributed(af::const_ref<double>(pixels.begin(), num_data))) {
          break;
        }
      }

      // Set rejected pixels as 'not background'
      for (std::size_t i = 0; i < num_data; ++i) {
        mask[indices[i]] |= shoebox::BackgroundUsed;
      }
    }

  private:
    std::size_t min_data_;
  };
}}}  // namespace dials::algorithms::background

#endif /* DIALS_ALGORITHMS_BACKGROUND_NORMAL_OUTLIER_REJECTOR_H */
