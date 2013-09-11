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

#include <dials/model/data/reflection.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using dials::model::Reflection;

  /**
   * Set the shoebox background value.
   * @param reflections The reflection list
   * @param value The value to set the background pixels to
   */
  void set_shoebox_background_value(
      af::ref<Reflection> reflections, double value) {
    for (std::size_t i = 0; i < reflections.size(); ++i) {
      af::ref< double, af::c_grid<3> > background =
        reflections[i].get_shoebox_background().ref();
      for (std::size_t j = 0; j < background.size(); ++j) {
        background[j] = value;
      }
    }
  }

}}

#endif /* DIALS_ALGORITHMS_BACKGROUND_DISCRIMINATOR_STRATEGY_H */
