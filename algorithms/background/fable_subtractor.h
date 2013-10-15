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

#include <omptbx/omp_or_stubs.h>
#include <dials/model/data/reflection.h>
#include <dials/error.h>
#include <dials/algorithms/background/poisson_discriminator.h>

namespace dials { namespace algorithms {

  using dials::model::Reflection;

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
    double operator()(const af::const_ref<double, af::c_grid<3> > &shoebox,
                      af::ref< int, af::c_grid<3> > mask) const {

      // Set which pixels belong in the background and which are spots
      discriminate_(shoebox.as_1d(), mask.as_1d());

      // Copy the background pixels into an array
      af::shared<double> pixels;
      for (std::size_t i = 0; i < mask.size(); ++i) {
        if (mask[i] & (shoebox::Background | shoebox::Valid)) {
          pixels.push_back(shoebox[i]);
        }
      }

      // Calculate the mean of the background pixels
      return mean(pixels.const_ref());
    }

    /**
     * Process the reflection
     * @params reflection The reflection
     */
    void operator()(Reflection &reflection) const {
      af::ref< int, af::c_grid<3> > mask = reflection.get_shoebox_mask().ref();
      af::ref< double, af::c_grid<3> > shoebox = reflection.get_shoebox().ref();
      af::ref< double, af::c_grid<3> > background =
        reflection.get_shoebox_background().ref();
      double value = this->operator()(shoebox, mask);
      for (std::size_t i = 0; i < background.size(); ++i) {
        background[i] = value;
      }
    }

    /**
     * Process the reflection list
     * @params reflections The list of reflections
     * @return Arrays of booleans True/False successful.
     */
    void operator()(af::ref<Reflection> reflections) const {
      #pragma omp parallel for
      for (std::size_t i = 0; i < reflections.size(); ++i) {
        try {
          if (reflections[i].is_valid()) {
            this->operator()(reflections[i]);
          }
        } catch(dials::error) {
          reflections[i].set_valid(false);
        }
      }
    }

    protected:

      PoissonDiscriminator discriminate_;
  };

}}

#endif /* DIALS_ALGORITHMS_BACKGROUND_FABLE_SUBTRACTOR_H */
