/*
 * poisson_discriminator.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_BACKGROUND_POISSON_DISCRIMINATOR_H
#define DIALS_ALGORITHMS_BACKGROUND_POISSON_DISCRIMINATOR_H

#include <cmath>
#include <scitbx/array_family/ref_reductions.h>
#include <scitbx/math/mean_and_variance.h>
#include <dials/array_family/sort_index.h>
#include <dials/model/data/reflection.h>
#include <dials/algorithms/shoebox/mask_code.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using std::sqrt;
  using std::pow;
  using scitbx::math::mean_and_variance;
  using dials::af::sort_index;
  using dials::model::Reflection;

  /**
   * Calculate the kth central moment.
   * @param data The data array
   * @param c The centre
   * @param k The number of the moment
   * @return The moment
   */
  inline
  double moment(const af::const_ref<double> &data, double c, std::size_t k) {
    std::size_t n = data.size();
    double m = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      m += pow((double)(data[i] - c), (int)k);
    }
    return m / n;
  }

  /**
   * Check if the data is poisson distributed.
   *
   * True is the absolute difference between the mean and sample variance
   * is less than a given number of standard deviations of the variance.
   *
   * @param data The array of pixel values
   * @param n_sigma The number of standard deviations
   * @returns True/False
   */
  inline
  bool is_poisson_distributed(const af::const_ref<double> &data, double n_sigma) {

    // Calculate the mean and standard deviation of the data
    mean_and_variance <double> mean_and_variance(data);
    double m1 = mean_and_variance.mean();
    double m2 = mean_and_variance.unweighted_sample_variance();

    // Estmate the variance of the variance and get the sdev
    double m4 = moment(data, m1, 4);

    // Estimate the standard deviation of the variance
    std::size_t n = data.size();
    double sdev = std::sqrt((m4 - m2 * m2 * (n - 3) / (n - 1)) / n);

    // Return True/False
    return std::abs(m2 - m1) <= n_sigma * sdev;
  }

  /**
   * A class that uses poisson distribution statistics to discriminate
   * between background and peak pixels in a reflection shoebox.
   */
  class PoissonDiscriminator {
  public:

    /** Initialise the class. */
    PoissonDiscriminator()
      : min_data_(10),
        n_sigma_(3.0) {}

    /** Initialise the class with parameters. */
    PoissonDiscriminator(std::size_t min_data, double n_sigma)
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
     * intensities are poisson distributed. If not then remove the pixel
     * with the highest intensity from the list and check again. Keep going
     * untill the list of pixels is poisson distributed, or the maximum
     * number of iterations is reached. The remaining pixels are classed
     * as background, the rest are peak.
     *
     * The test performed is to check that the absolute difference between
     * the sample variance and the mean is less than a certain number of
     * standard deviations of the variance.
     *
     * @params shoebox The shoebox profile
     * @params mask The shoebox mask
     */
    void operator()(const af::const_ref<double> &shoebox, af::ref<int> mask) const {

      // Ensure data is correctly sized.
      DIALS_ASSERT(shoebox.size() == mask.size());

      // Copy valid pixels and indices into list
      af::shared<int> indices;
      for (std::size_t i = 0; i < shoebox.size(); ++i) {
        if (mask[i] & shoebox::Valid) {
          indices.push_back(i);
        }
      }

      // Check we have enough data
      DIALS_ASSERT(indices.size() >= min_data_);

      // Check we have enough data
      DIALS_ASSERT(indices.size() >= min_data_);

      // Sort the pixels into ascending intensity order
      sort_index(indices.begin(), indices.end(), shoebox.begin());
      af::shared<double> pixels(indices.size());
      for (std::size_t i = 0; i < indices.size(); ++i) {
        pixels[i] = (double)shoebox[indices[i]];
      }

      // Check if the data is poissson distributed. If it is not, then remove
      // a value of high intensity and keep looping until it is. If the number
      // of iterations exceeds the maximum then exit the loop.
      std::size_t num_data = pixels.size();
      for (; num_data > min_data_; --num_data) {
        if (is_poisson_distributed(af::const_ref<double>(
            pixels.begin(), num_data), n_sigma_)) {
          break;
        }
      }

      // Set all the rejected pixels as peak pixels and all the accepted
      // pixels as background pixels
      for (std::size_t i = 0; i < num_data; ++i) {
        mask[indices[i]] |= shoebox::Background;
      }
      for (std::size_t i = num_data; i < indices.size(); ++i) {
        mask[indices[i]] |= shoebox::Foreground;
      }
    }

    /**
     * Process just a shoebox and return a mask
     * @param shoebox The shoebox profile
     * @return The mask
     */
    af::shared<int> operator()(const af::const_ref<double> &shoebox) const {
      af::shared<int> mask(shoebox.size(), shoebox::Valid);
      af::ref<int> mask_ref = mask.ref();
      this->operator()(shoebox, mask_ref);
      return mask;
    }

    /**
     * Process the reflection
     * @param reflection The reflection
     */
    void operator()(Reflection &reflection) const {
      af::const_ref<double> shoebox = reflection.get_shoebox().const_ref().as_1d();
      af::ref<int> mask = reflection.get_shoebox_mask().ref().as_1d();
      this->operator()(shoebox, mask);
    }

  private:

    std::size_t min_data_;
    double n_sigma_;
  };
}}

#endif /* DIALS_ALGORITHMS_BACKGROUND_POISSON_DISCRIMINATOR_H */
