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

#include <dials/model/data/shoebox.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using dials::model::Shoebox;

  /**
   * Set the shoebox background value.
   * @param reflections The reflection list
   * @param value The value to set the background pixels to
   */
  template <typename FloatType>
  void set_shoebox_background_value(af::ref<Shoebox<FloatType> > shoeboxes,
                                    FloatType value) {
    for (std::size_t i = 0; i < shoeboxes.size(); ++i) {
      af::ref<FloatType, af::c_grid<3> > background = shoeboxes[i].background.ref();
      for (std::size_t j = 0; j < background.size(); ++j) {
        background[j] = value;
      }
    }
  }

}}  // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_BACKGROUND_DISCRIMINATOR_STRATEGY_H */
