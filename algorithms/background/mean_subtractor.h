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

#include <dials/algorithms/background/subtractor_strategy.h>

namespace dials { namespace algorithms {

  using scitbx::af::flex_bool;
  using dials::model::ReflectionList;

  /** Class for background subtraction by mean background value */
  class MeanSubtractor : public SubtractorStrategy {
  public:

    /** Initialise the class. */
    MeanSubtractor() {}

    /**
     * Process the reflection list
     * @params reflections The list of reflections
     * @return Arrays of booleans True/False successful.
     */
    virtual flex_bool operator()(ReflectionList &reflections) const {
      return flex_bool();
    }
  };

}}

#endif /* DIALS_ALGORITHMS_BACKGROUND_MEAN_SUBTRACTOR_H */
