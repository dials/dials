/*
 * subtract_background.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_INTEGRATION_SUBTRACT_BACKGROUND_H
#define DIALS_ALGORITHMS_INTEGRATION_SUBTRACT_BACKGROUND_H

#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/flex_types.h>
#include <dials/model/data/reflection.h>
#include <dials/error.h>
#include "background_intensity.h"

namespace dials { namespace algorithms {

  using scitbx::af::flex_bool;
  using scitbx::af::flex_int;
  using scitbx::af::flex_double;
  using dials::model::Reflection;
  using dials::model::ReflectionList;

  /** A class to subtract the background intensity from the reflection profile */
  class SubtractBackground {

  public:

    /**
     * Initialise the class with parameters
     * @param delta The deviation from normal
     * @param max_iter The maximum numner of iterations as a fraction of the
     *        elements in the input data array
     * @param min_pixels The minimum number of pixels needed to calculate
     *        the background intensity
     */
    SubtractBackground(int min_pixels = 10, double n_sigma = -1)
      : min_pixels_(min_pixels),
        n_sigma_(n_sigma) {}

    /**
     * Calculate the background intensity for a single reflection and subtract
     * it from the image pixel values.
     *
     * @todo In the XDS paper, the background intensity value is over estimated
     *       for strong reflections and is adjusted using the modelled
     *       intensity profile in the xds frame. This needs to be done.
     *
     * @param reflection The reflection
     */
    void operator()(Reflection &reflection)
    {
      // Get the shoebox
      flex_int shoebox = reflection.get_shoebox();
      flex_int mask    = reflection.get_shoebox_mask();

      // Ensure data is correctly sized.
      DIALS_ASSERT(shoebox.size() == mask.size());

      // Copy valid pixels into list
      std::size_t count = 0;
      flex_double shoebox_pixels = flex_double(shoebox.size());
      for (std::size_t i = 0; i < shoebox.size(); ++i) {
        if (mask[i]) {
          shoebox_pixels[count++] = (double)shoebox[i];
        }
      }
      shoebox_pixels.resize(count);

      // Estimate the background intensity
      double background = background_intensity(shoebox_pixels,
          min_pixels_, n_sigma_);

      // Subtract background from shoebox profile
      for (std::size_t i = 0; i < shoebox.size(); ++i) {
        shoebox[i] -= (int)background;
      }
    }

    /**
     * Subtract the background for all reflections
     * @param reflections The array of reflections
     * @returns The a boolean array containing the status for each reflection.
     *          True/False was the background successfully subtracted
     */
    flex_bool operator()(ReflectionList &reflections) {
      flex_bool result(reflections.size());
      for (int i = 0; i < reflections.size(); ++i) {
        try {
          this->operator()(reflections[i]);
          result[i] = true;
        } catch(error) {
          result[i] = false;
        }
      }
      return result;
    }

  private:

    int min_pixels_;
    double n_sigma_;
  };

}} // namespace = dials::algorithms

#endif // DIALS_ALGORITHMS_INTEGRATION_SUBTRACT_BACKGROUND_H
