/*
 * reflection_profile_helpers.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_INTEGRATION_REFLECTION_PROFILE_HELPERS_H
#define DIALS_ALGORITHMS_INTEGRATION_REFLECTION_PROFILE_HELPERS_H

#include <scitbx/array_family/flex_types.h>
#include <dials/model/data/reflection.h>
#include <dials/algorithms/shoebox/mask_code.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using scitbx::af::int2;
  using scitbx::af::int6;
  using scitbx::af::flex_bool;
  using scitbx::af::flex_int;
  using scitbx::af::flex_double;
  using scitbx::af::flex_grid;
  using dials::model::Reflection;
  using dials::model::ReflectionList;

  /**
   * Copy the pixels from a single image into the profile arrays of
   * the reflections that are recorded on this image.
   *
   * @param image The 2D image
   * @param array_index The array index of the image
   * @param index The reflection indices
   * @param reflections The reflection list.
   * @param gain_map The detector gain map
   * @param dark_map The detector dark map
   */
  inline
  void assign_strong_spots(const flex_bool &mask, int array_index,
      const flex_int &index, ReflectionList &reflections) {
    for (std::size_t i = 0; i < index.size(); ++i) {

      // Get a reference to a reflection
      Reflection &r = reflections[index[i]];
      if (r.is_valid()) {
        int6 bounding_box = r.get_bounding_box();
        int i0 = bounding_box[0], i1 = bounding_box[1];
        int j0 = bounding_box[2], j1 = bounding_box[3];
        int k0 = bounding_box[4];

        int k = array_index - k0;

        // Get the image size
        flex_int::index_type mask_size = mask.accessor().all();

        // Readjust the area to loop over to ensure we're within image bounds
        int jj0 = j0 >= 0 ? j0 : 0;
        int ii0 = i0 >= 0 ? i0 : 0;
        int jj1 = j1 <= mask_size[0] ? j1 : mask_size[0];
        int ii1 = i1 <= mask_size[1] ? i1 : mask_size[1];

        // Get the reflection profile
        flex_int shoebox_mask = r.get_shoebox_mask();

        // Assign strong pixels (if any)
        std::size_t strong_count = 0;
        for (int jj = jj0; jj < jj1; ++jj) {
          for (int ii = ii0; ii < ii1; ++ii) {
            int j = jj - j0;
            int i = ii - i0;
            if (shoebox_mask(k, j, i) & shoebox::Valid) {
              bool strong = mask(jj, ii);
              if (strong) {
                shoebox_mask(k, j, i) |= shoebox::Strong;
                strong_count++;
              }
            }
          }
        }

        // If strong count is > 0 then spot is strong
        if (strong_count > 0) {
          r.set_strong(true);
        }
      }
    }
  }


}} // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_INTEGRATION_REFLECTION_PROFILE_HELPERS_H */
