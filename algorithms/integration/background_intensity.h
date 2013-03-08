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
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/ref_reductions.h>
#include <boost/math/special_functions/erf.hpp>
#include <scitbx/math/mean_and_variance.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using boost::math::erfc;
  using scitbx::af::shared;
  using scitbx::af::const_ref;
  using scitbx::af::ref;
  using scitbx::af::min;
  using scitbx::af::max;
  using scitbx::af::mean;
  using scitbx::math::mean_and_variance;

  /**
   * Get the expected number of standard deviations based on the number of
   * observations. Given by erf(x / sqrt(2)) = 1 - 1 / N
   * @param n_obs The number of observations
   * @returns The expected number of standard deviations
   */
  inline
  double expected_n_sigma(int n_obs) {
    return sqrt(2.0) * erfc(1.0 - (1.0 / n_obs));
  }

  /**
   * Get the maximum number of standard deviations in the range of data
   * @param n_obs The number of observations
   * @returns The expected number of standard deviations
   */
  inline
  double maximum_n_sigma(const const_ref<double> &data) {

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
    double max_n_sigma = maximum_n_sigma(data);

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
    double n_sigma = expected_n_sigma(data.size());

    // Check if data is normally distributed using sigma value
    return is_normally_distributed(data, n_sigma);
  }

  /**
   * Functor to compare in sort_index.
   */
  template <class T>
  struct index_less {
    index_less(const T &v) : v_(v) {}

    template <class IndexType>
    bool operator() (const IndexType& x, const IndexType& y) const {
      return v_[x] < v_[y];
    }
    const T &v_;
  };

  /**
   * Given a vector return a sorted list of indices.
   * @param v The list of values
   * @returns A sorted list of indices
   */
  template <typename T>
  shared<std::size_t> sort_index(const const_ref<T> &v) {

    // initialize original index locations
    shared<std::size_t> index(v.size());
    for (size_t i = 0; i != index.size(); ++i) {
      index[i] = i;
    }

    // sort indexes based on comparing values in v
    std::sort(index.begin(), index.end(), index_less<const_ref<T> >(v));

    // Return indices
    return index;
  }

  /*
   * Calculate the pixels contributing to the background intensity.
   *
   * Sort the pixels in order of ascending intensity. Then check if the
   * intensities are normally distributed. If not then remove the pixel
   * with the highest intensity from the list and check again. Keep going
   * untill the list of pixels is normally distributed, or the maximum
   * number of iterations is reached.
   *
   * @param data The list of pixels
   * @param min_data The minimum number of pixels needed
   * @param n_sigma The number of standard deviations to consider normal
   * @returns The list of pixels contributing to the background
   */
  shared<std::size_t> background_pixels(const const_ref<double> &data,
      int min_data, double n_sigma) {

    // Check we have enough data
    DIALS_ASSERT(min_data > 0);
    DIALS_ASSERT(data.size() >= min_data);

    // Sort the data and return sorted indices
    shared<std::size_t> index = sort_index(data);
    shared<double> sorted_data(data.size());
    for (std::size_t i = 0; i < index.size(); ++i) {
      sorted_data[i] = data[index[i]];
    }

    // Check if the data is normally distributed. If it is not, then remove
    // a value of high intensity and keep looping until it is. If the number
    // of iterations exceeds the maximum then exit the loop.
    std::size_t num_data = data.size();
    for (; num_data > min_data; --num_data) {
      if (is_normally_distributed(const_ref<double>(
          sorted_data.begin(), num_data), n_sigma)) {
        break;
      }
    }

    // Resize index to first n elements
    index.resize(num_data);

    // Return the indices
    return index;
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
   * @param pixels The list of pixels
   * @returns The background intensity value
   */
  inline
  double background_intensity(const shared<double> &data, int min_data,
      double n_sigma)
  {
    // Check we have enough data
    DIALS_ASSERT(data.size() >= min_data);

    shared<double> sorted_data(data);

    // Sort the pixels into ascending intensity order
    std::sort(sorted_data.begin(), sorted_data.end());

    // Check if the data is normally distributed. If it is not, then remove
    // a value of high intensity and keep looping until it is. If the number
    // of iterations exceeds the maximum then exit the loop.
    std::size_t num_data = data.size();
    for (; num_data > min_data; --num_data) {
      if (is_normally_distributed(const_ref<double>(
          sorted_data.begin(), num_data), n_sigma)) {
        break;
      }
    }

    // Return the mean of the remaining pixels as the background intensity
    return mean(const_ref<double>(sorted_data.begin(), num_data));
  }
}}

#endif /* DIALS_ALGORITHMS_INTEGRATION_BACKGROUND_INTENSITY_H */
