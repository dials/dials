/*
 * subtractor_strategy.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_BACKGROUND_SUBTRACTOR_STRATEGY_H
#define DIALS_ALGORITHMS_BACKGROUND_SUBTRACTOR_STRATEGY_H

#include <scitbx/array_family/flex_types.h>
#include <dials/model/data/reflection.h>

namespace dials { namespace algorithms {

  using scitbx::af::flex_bool;
  using dials::model::ReflectionList;

  /** Base class for background subtraction strategies */
  class SubtractorStrategy {
  public:

    /**
     * Process the reflection list
     * @params reflections The list of reflections
     * @return Arrays of booleans True/False successful.
     */
    virtual flex_bool operator()(ReflectionList &reflections) const = 0;
  };

}}

#endif /* DIALS_ALGORITHMS_BACKGROUND_SUBTRACTOR_STRATEGY_H */
