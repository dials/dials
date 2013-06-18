/*
 * helpers.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_BACKGROUND_HELPERS_H
#define DIALS_ALGORITHMS_BACKGROUND_HELPERS_H

#include <scitbx/array_family/flex_types.h>
#include <dials/model/data/reflection.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using scitbx::af::flex_double;
  using dials::model::Reflection;
  using dials::model::ReflectionList;

  /**
   * Set the shoebox background value.
   * @param reflections The reflection list
   * @param value The value to set the background pixels to
   */
  void set_shoebox_background_value(ReflectionList &reflections, double value) {
    for (std::size_t i = 0; i < reflections.size(); ++i) {
      flex_double background = reflections[i].get_shoebox_background();
      for (std::size_t j = 0; j < background.size(); ++j) {
        background[j] = value;
      }
    }
  }

}}

#endif /* DIALS_ALGORITHMS_BACKGROUND_DISCRIMINATOR_STRATEGY_H */
