/*
 * index_of_dispersion_discriminator.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_BACKGROUND_INDEX_OF_DISPERSION_DISCRIMINATOR_H
#define DIALS_ALGORITHMS_BACKGROUND_INDEX_OF_DISPERSION_DISCRIMINATOR_H

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
   * Check if the data is poisson distributed. If the index of dispersion
   * (which should be 1) is greater than n_sigma * the standard deviation of
   * the expected IOD, then data is not poisson distributed.
   * @param data The data to use
   * @param n_sigma The number of standard deviations.
   * @returns True/False the data is poisson distributed.
   */
  inline
  bool is_poisson_distributed_iod(const const_ref<double> &data, double n_sigma)
  {
    // Calculate the mean and variance of the data
    mean_and_variance <double> mean_and_variance(data);
    double mean = mean_and_variance.mean();
    double var  = mean_and_variance.unweighted_sample_variance();
    double count = mean_and_variance.sum_weights();

    // Check the variance
    if (var <= 0) {
      return true;
    }

    // Calculate the index of dispersion and check if it is within the bounds
    return (var / mean) <= 1.0 + n_sigma * sqrt(2.0 / (count - 1));
  }

  /**
   * A class that uses normal distribution statistics to discriminate
   * between background and peak pixels in a reflection shoebox.
   */
  class IndexOfDispersionDiscriminator : public DiscriminatorStrategy {
  public:

    /** Initialise the class. */
    IndexOfDispersionDiscriminator()
      : min_data_(10),
        n_sigma_(3.0) {}

    /** Initialise the class with parameters. */
    IndexOfDispersionDiscriminator(std::size_t min_data, double n_sigma)
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
     * untill the list of pixels is poisson distributed, or the maximum
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
        if (is_poisson_distributed_iod(const_ref<double>(
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

#endif /* DIALS_ALGORITHMS_BACKGROUND_INDEX_OF_DISPERSION_DISCRIMINATOR_H */
