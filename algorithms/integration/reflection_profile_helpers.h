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
   * Allocate the memory for reflection profiles (on the detector side)
   * and the reflection mask.
   * @param reflections The reflection list.
   */
  inline
  ReflectionList allocate_reflection_profiles(ReflectionList &reflections,
      double shoebox_default = 0, int shoebox_mask_default = 1,
      double shoebox_background_default = 0) {
    for (std::size_t i = 0; i < reflections.size(); ++i) {
      Reflection &r = reflections[i];
      if (r.is_valid()) {
        int size_z = r.get_bounding_box()[5] - r.get_bounding_box()[4];
        int size_y = r.get_bounding_box()[3] - r.get_bounding_box()[2];
        int size_x = r.get_bounding_box()[1] - r.get_bounding_box()[0];
        DIALS_ASSERT(size_z > 0 && size_y > 0 && size_x > 0);
        r.set_shoebox(flex_double(flex_grid<>(size_z, size_y, size_x),
          shoebox_default));
        r.set_shoebox_mask(flex_int(flex_grid<>(size_z, size_y, size_x),
          shoebox_mask_default));
        r.set_shoebox_background(flex_double(flex_grid<>(size_z, size_y, size_x),
          shoebox_background_default));
      }
    }
    return reflections;
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
  void copy_single_image_pixels(const flex_int &image, int array_index,
      const flex_int &index, ReflectionList &reflections,
      flex_double gain_map, flex_int dark_map) {
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
        flex_int::index_type image_size = image.accessor().all();

        // Readjust the area to loop over to ensure we're within image bounds
        int jj0 = j0 >= 0 ? j0 : 0;
        int ii0 = i0 >= 0 ? i0 : 0;
        int jj1 = j1 <= image_size[0] ? j1 : image_size[0];
        int ii1 = i1 <= image_size[1] ? i1 : image_size[1];

        // Get the reflection profile
        flex_double profile = r.get_shoebox();

        // Copy the image pixels
        for (int jj = jj0; jj < jj1; ++jj) {
          for (int ii = ii0; ii < ii1; ++ii) {
            int j = jj - j0;
            int i = ii - i0;
            profile(k, j, i) = gain_map(jj, ii) * (
                image(jj, ii) - dark_map(jj, ii));
          }
        }

        // Set the reflection profile
        r.set_shoebox(profile);
      }
    }
  }

  /**
   * Construct a mask of pixels to use in the determination of strong spots.
   * Pixels that are to be used are returned as 1, otherwise as 0.
   * @param mask The mask input
   * @param array_index The index of the image
   * @param index The indices of reflections on this image
   * @param reflections The reflection array
   * @param kernel_size The size of the kernel to use around the shoebox
   */
  inline
  flex_bool construct_image_mask_from_shoeboxes(const flex_bool &mask,
      int array_index, const flex_int &index, const ReflectionList &reflections,
      int2 kernel_size) {

    // Create the resulting mask
    flex_bool result(mask.accessor(), false);

    // Set all the shoebox pixels
    for (std::size_t i = 0; i < index.size(); ++i) {

      // Get a reference to a reflection
      const Reflection &r = reflections[index[i]];
      if (r.is_valid()) {
        int6 bounding_box = r.get_bounding_box();
        int i0 = bounding_box[0] - kernel_size[0];
        int i1 = bounding_box[1] + kernel_size[0];
        int j0 = bounding_box[2] - kernel_size[1];
        int j1 = bounding_box[3] + kernel_size[1];

        // Get the image size
        flex_int::index_type mask_size = mask.accessor().all();

        // Readjust the area to loop over to ensure we're within image bounds
        int jj0 = j0 >= 0 ? j0 : 0;
        int ii0 = i0 >= 0 ? i0 : 0;
        int jj1 = j1 <= mask_size[0] ? j1 : mask_size[0];
        int ii1 = i1 <= mask_size[1] ? i1 : mask_size[1];

        // Copy the image pixels
        for (int jj = jj0; jj < jj1; ++jj) {
          for (int ii = ii0; ii < ii1; ++ii) {
            result(jj, ii) = true;
          }
        }
      }
    }

    // Mask the pixels
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = result[i] && mask[i];
    }

    // Return the resulting mask
    return result;
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
            if (shoebox_mask(k, j, i) != 0) {
              bool strong = mask(jj, ii);
              if (strong) {
                shoebox_mask(k, j, i) |= 2;
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
