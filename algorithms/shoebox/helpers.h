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
#ifndef DIALS_ALGORITHMS_SHOEBOX_HELPERS_H
#define DIALS_ALGORITHMS_SHOEBOX_HELPERS_H

#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/tiny_types.h>
#include <dials/model/data/reflection.h>
#include <dials/algorithms/shoebox/mask_code.h>
#include <dials/error.h>

namespace dials { namespace algorithms { namespace shoebox {

  using scitbx::af::flex_grid;
  using scitbx::af::flex_int;
  using scitbx::af::flex_double;
  using scitbx::af::int6;
  using dials::model::Reflection;
  using dials::model::ReflectionList;

  /**
   * Allocate the profiles in the reflection list
   * @param reflections The reflection list
   */
  inline
  void allocate(ReflectionList &reflections) {
    // Allocate all the reflection profiles
    for (std::size_t i = 0; i < reflections.size(); ++i) {
      Reflection &r = reflections[i];
      if (r.is_valid()) {
        int6 bbox = r.get_bounding_box();
        int size_z = bbox[5] - bbox[4];
        int size_y = bbox[3] - bbox[2];
        int size_x = bbox[1] - bbox[0];
        DIALS_ASSERT(size_z > 0 && size_y > 0 && size_x > 0);
        flex_grid<> accessor(size_z, size_y, size_x);
        r.set_shoebox(flex_double(accessor, 0.0));
        r.set_shoebox_mask(flex_int(accessor, shoebox::Valid));
        r.set_shoebox_background(flex_double(accessor, 0.0));
      }
    }
  }

  /**
   * Deallocate the profiles in the reflection list
   * @param reflections The reflection list
   */
  inline
  void deallocate(ReflectionList &reflections) {
    // Delete all the reflection profiles from memory
    for (std::size_t i = 0; i < reflections.size(); ++i) {
      reflections[i].set_shoebox(flex_double());
      reflections[i].set_shoebox_mask(flex_int());
      reflections[i].set_shoebox_background(flex_double());
      reflections[i].set_transformed_shoebox(flex_double());
    }
  }

}}} // namespace dials::algorithms::shoebox

#endif /* DIALS_ALGORITHMS_SHOEBOX_HELPERS_H */
