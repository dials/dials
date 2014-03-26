/*
 * xds_subtractor.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_BACKGROUND_XDS_SUBTRACTOR_H
#define DIALS_ALGORITHMS_BACKGROUND_XDS_SUBTRACTOR_H

#include <omptbx/omp_or_stubs.h>
#include <dials/error.h>
#include <dials/model/data/shoebox.h>
#include <dials/algorithms/background/normal_discriminator.h>

namespace dials { namespace algorithms {

  using model::Shoebox;

  /** The XDS background subtraction algoritm */
  class XdsSubtractor {
  public:

    /**
     * Initialise the algorithm.
     * @param min_data The minimum number of pixels to use
     * @param n_sigma The number of standard deviations
     */
    XdsSubtractor(std::size_t min_data):
      discriminate_(min_data) {}

    /**
     * Process the shoebox
     * @params shoebox The shoebox
     * @params mask The shoebox mask
     * @returns The background value
     */
    template <typename FloatType>
    FloatType operator()(
        const af::const_ref< FloatType, af::c_grid<3> > &shoebox,
        af::ref< int, af::c_grid<3> > mask) const {

      // Set which pixels belong in the background and which are spots
      discriminate_(shoebox.as_1d(), mask.as_1d());

      // Copy the background pixels into an array
      af::shared<FloatType> pixels;
      for (std::size_t i = 0; i < mask.size(); ++i) {
        if (mask[i] & shoebox::Valid && mask[i] & shoebox::BackgroundUsed) {
          pixels.push_back(shoebox[i]);
        }
      }

      // Calculate the mean of the background pixels
      return mean(pixels.const_ref());
    }

  protected:

    NormalDiscriminator discriminate_;
  };

}}

#endif /* DIALS_ALGORITHMS_BACKGROUND_XDS_SUBTRACTOR_H */
