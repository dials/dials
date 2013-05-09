/*
 * discriminator_strategy.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_BACKGROUND_DISCRIMINATOR_STRATEGY_H
#define DIALS_ALGORITHMS_BACKGROUND_DISCRIMINATOR_STRATEGY_H

#include <scitbx/array_family/flex_types.h>
#include <dials/model/data/reflection.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using scitbx::af::flex_bool;
  using dials::model::Reflection;
  using dials::model::ReflectionList;

  /** Base class for pixel discrimination strategies */
  class DiscriminatorStrategy {
  public:

    /**
     * Process the reflection
     * @params reflection The reflection
     */
    virtual void operator()(Reflection &reflection) const = 0;

    /**
     * Process the reflection list
     * @params reflections The list of reflections
     * @return Arrays of booleans True/False successful.
     */
    flex_bool operator()(ReflectionList &reflections) const {
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
  };

}}

#endif /* DIALS_ALGORITHMS_BACKGROUND_DISCRIMINATOR_STRATEGY_H */
