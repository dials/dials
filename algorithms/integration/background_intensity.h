/*
 * background_intensity.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_INTEGRATION_BACKGROUND_INTENSITY_H
#define DIALS_ALGORITHMS_INTEGRATION_BACKGROUND_INTENSITY_H

#include <algorithm>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/ref_reductions.h>
#include <boost/math/special_functions/erf.hpp>
#include <scitbx/math/mean_and_variance.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using boost::math::erfc;
  using scitbx::af::const_ref;
  using scitbx::af::ref;
  using scitbx::af::min;
  using scitbx::af::max;
  using scitbx::af::mean;
  using scitbx::math::mean_and_variance;
  
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
  bool is_normally_distributed(const const_ref<double> &data, int n_sigma) {
 
    // Calculate the mean and standard deviation of the data
    mean_and_variance <double> mean_and_variance(data);
    double mean = mean_and_variance.mean();
    double sdev = mean_and_variance.unweighted_sample_standard_deviation();
    DIALS_ASSERT(sdev > 0);

    // Calculate the min/max of the data
    double mind = min(data);
    double maxd = max(data);
    
    // Calculate t-statistic of min/max
    double min_n_sigma = (mean - mind) / sdev;
    double max_n_sigma = (maxd - mean) / sdev;
    
    // return whether within required sigma
    return max_n_sigma < n_sigma && min_n_sigma < n_sigma;
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
    double n_sigma = sqrt(2.0) * erfc(1.0 - (1.0 / data.size()));
 
    // Check if data is normally distributed using sigma value
    return is_normally_distributed(data, n_sigma);
  }

  /*
   * Calculate the background intensity.
   *
   * Sort the pixels in order of ascending intensity. Then check if the
   * intensities are normally distributed. If not then remove the pixel
   * with the highest intensity from the list and check again. Keep going
   * untill the list of pixels is normally distributed, or the maximum
   * number of iterations is reached. Return the mean of the values as the
   * background intensity.
   *
   * This function modifies the input data
   *
   * @param pixels The list of pixels
   * @returns The background intensity value
   */
  inline
  double background_intensity(ref<double> &data, int min_data, double n_sigma)
  {
    // Check we have enough data
    DIALS_ASSERT(data.size() >= min_data);

    // Sort the pixels into ascending intensity order
    std::sort(data.begin(), data.end());

    // Check if the data is normally distributed. If it is not, then remove
    // a value of high intensity and keep looping until it is. If the number
    // of iterations exceeds the maximum then exit the loop.
    int num_data = data.size();
    for (; num_data >= min_data; --num_data) {
      if (is_normally_distributed(const_ref<double>(
          data.begin(), num_data), n_sigma)) {
        break;
      }
    }

    // Return the mean of the remaining pixels as the background intensity
    return mean(const_ref<double>(data.begin(), num_data));
  }
}}

#endif /* DIALS_ALGORITHMS_INTEGRATION_BACKGROUND_INTENSITY_H */

