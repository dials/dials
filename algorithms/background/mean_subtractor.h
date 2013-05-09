/*
 * mean_subtractor.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_BACKGROUND_MEAN_SUBTRACTOR_H
#define DIALS_ALGORITHMS_BACKGROUND_MEAN_SUBTRACTOR_H

#include <scitbx/array_family/ref_reductions.h>
#include <dials/algorithms/background/subtractor_strategy.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using scitbx::af::shared;
  using scitbx::af::flex_int;
  using scitbx::af::flex_double;
  using scitbx::af::mean;

  /** Class for background subtraction by mean background value */
  class MeanSubtractor : public SubtractorStrategy {
  public:

    /** Initialise the class. */
    MeanSubtractor() {}

    /**
     * Process the reflection list
     * @params reflection The reflections
     */
    virtual void operator()(Reflection &reflection) const {

      // Get the shoebox
      flex_int shoebox = reflection.get_shoebox();

      // Calculate the background intensity
      int background = background_intensity(reflection);

      // Subtract background from shoebox profile
      for (std::size_t i = 0; i < shoebox.size(); ++i) {
        shoebox[i] -= background;
      }
    }

  protected:

    /**
     * Calculate the background intensity.
     * @param reflection The reflection
     * @return The background intensity.
     */
    virtual int background_intensity(Reflection &reflection) const {

      // Get the shoebox and mask and ensure the correct size
      flex_int mask = reflection.get_shoebox_mask();
      flex_int shoebox = reflection.get_shoebox();
      DIALS_ASSERT(mask.size() == shoebox.size());

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
 };

}}

#endif /* DIALS_ALGORITHMS_BACKGROUND_MEAN_SUBTRACTOR_H */
