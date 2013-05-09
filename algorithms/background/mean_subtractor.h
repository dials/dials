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
    MeanSubtractor(std::size_t min_pixels)
      : min_pixels_(min_pixels) {
      DIALS_ASSERT(min_pixels_ > 0);
    }

    /**
     * Process the reflection list
     * @params reflection The reflections
     */
    virtual void operator()(Reflection &reflection) const {

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

      // Ensure that we have enough pixels
      DIALS_ASSERT(pixels.size() > min_pixels_);

      // Calculate the mean of the background pixels
      int background = (int)mean(pixels.const_ref());

      // Subtract background from shoebox profile
      for (std::size_t i = 0; i < shoebox.size(); ++i) {
        shoebox[i] -= background;
      }
    }

  private:
    std::size_t min_pixels_;
  };

}}

#endif /* DIALS_ALGORITHMS_BACKGROUND_MEAN_SUBTRACTOR_H */
