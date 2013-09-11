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

#include <dials/model/data/reflection.h>
#include <dials/algorithms/shoebox/mask_code.h>
#include <dials/error.h>

namespace dials { namespace algorithms { namespace shoebox {

  using scitbx::af::int2;
  using scitbx::af::int6;
  using dials::model::Reflection;

  /**
   * Allocate the profiles in the reflection list
   * @param reflections The reflection list
   * @param mask_default The default mask value
   */
  inline
  void allocate(af::ref<Reflection> reflections, int mask_default) {
    // Allocate all the reflection profiles
    for (std::size_t i = 0; i < reflections.size(); ++i) {
      Reflection &r = reflections[i];
      if (r.is_valid()) {
        int6 bbox = r.get_bounding_box();
        int size_z = bbox[5] - bbox[4];
        int size_y = bbox[3] - bbox[2];
        int size_x = bbox[1] - bbox[0];
        DIALS_ASSERT(size_z > 0 && size_y > 0 && size_x > 0);
        af::c_grid<3> accessor(size_z, size_y, size_x);
        r.set_shoebox(af::versa< double, af::c_grid<3> >(accessor, 0.0));
        r.set_shoebox_mask(af::versa< int, af::c_grid<3> >(accessor, mask_default));
        r.set_shoebox_background(af::versa< double, af::c_grid<3> >(accessor, 0.0));
      }
    }
  }

  /**
   * Allocate the profiles in the reflection list
   * @param reflections The reflection list
   */
  inline
  void allocate(af::ref<Reflection> reflections) {
    allocate(reflections, shoebox::Valid);
  }

  /**
   * Deallocate the profiles in the reflection list
   * @param reflections The reflection list
   */
  inline
  void deallocate(af::ref<Reflection> reflections) {
    // Delete all the reflection profiles from memory
    for (std::size_t i = 0; i < reflections.size(); ++i) {
      reflections[i].set_shoebox(af::versa< double, af::c_grid<3> >());
      reflections[i].set_shoebox_mask(af::versa< int, af::c_grid<3> >());
      reflections[i].set_shoebox_background(af::versa< double, af::c_grid<3> >());
      reflections[i].set_transformed_shoebox(af::versa< double, af::c_grid<3> >());
    }
  }

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
  void assign_strong_spots(const af::const_ref< bool, af::c_grid<2> > &mask,
      int array_index, const af::const_ref<int> &index,
      af::ref<Reflection> reflections) {
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
        int2 mask_size = mask.accessor();

        // Readjust the area to loop over to ensure we're within image bounds
        int jj0 = j0 >= 0 ? j0 : 0;
        int ii0 = i0 >= 0 ? i0 : 0;
        int jj1 = j1 <= mask_size[0] ? j1 : mask_size[0];
        int ii1 = i1 <= mask_size[1] ? i1 : mask_size[1];

        // Get the reflection profile
        af::ref< int, af::c_grid<3> > shoebox_mask = r.get_shoebox_mask().ref();

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

}}} // namespace dials::algorithms::shoebox

#endif /* DIALS_ALGORITHMS_SHOEBOX_HELPERS_H */
