/*
 * fable_subtractor.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_BACKGROUND_FABLE_SUBTRACTOR_H
#define DIALS_ALGORITHMS_BACKGROUND_FABLE_SUBTRACTOR_H

#include <scitbx/array_family/flex_types.h>
#include <dials/model/data/reflection.h>
#include <dials/error.h>
#include <dials/algorithms/background/poisson_discriminator.h>

namespace dials { namespace algorithms {

  using scitbx::af::flex_bool;
  using dials::model::Reflection;
  using dials::model::ReflectionList;

  /** The fable background subtraction algorithm */
  class FableSubtractor {
  public:

    /**
     * Initialise the algorithm.
     * @param min_data The minimum number of pixels to use
     * @param n_sigma The number of standard deviations
     */
    FableSubtractor(std::size_t min_data, double n_sigma):
      discriminate_(min_data, n_sigma) {}

    /**
     * Process the shoebox
     * @params shoebox The shoebox
     * @params mask The shoebox mask
     * @returns The background value
     */
    int operator()(const flex_int &shoebox, flex_int &mask) const {

      // Set which pixels belong in the background and which are spots
      discriminate_(shoebox, mask);

      // Copy the background pixels into an array
      shared<double> pixels;
      for (std::size_t i = 0; i < mask.size(); ++i) {
        if (mask[i] & (1 << 0)) {
          pixels.push_back(shoebox[i]);
        }
      }

      // Calculate the mean of the background pixels
      return (int)mean(pixels.const_ref());
    }

    /**
     * Process the reflection
     * @params reflection The reflection
     */
    void operator()(Reflection &reflection) const {
      flex_int background = reflection.get_shoebox_background();
      int value = this->operator()(
          reflection.get_shoebox(),
          reflection.get_shoebox_mask());
      for (std::size_t i = 0; i < background.size(); ++i) {
        background[i] = value;
      }
    }

    /**
     * Process the reflection list
     * @params reflections The list of reflections
     * @return Arrays of booleans True/False successful.
     */
    void operator()(ReflectionList &reflections) const {
      for (int i = 0; i < reflections.size(); ++i) {
        try {
          this->operator()(reflections[i]);
        } catch(dials::error) {
          continue;
        }
      }
    }

    protected:

      PoissonDiscriminator discriminate_;
  };

}}

#endif /* DIALS_ALGORITHMS_BACKGROUND_FABLE_SUBTRACTOR_H */
