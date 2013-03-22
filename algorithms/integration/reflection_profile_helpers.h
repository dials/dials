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

namespace dials { namespace algorithms {

  using scitbx::af::int6;
  using scitbx::af::flex_int;
  using scitbx::af::flex_grid;
  using dials::model::Reflection;
  using dials::model::ReflectionList;

  /**
   * Allocate the memory for reflection profiles (on the detector side).
   * @param reflections The reflection list.
   */
  inline
  void allocate_reflection_profiles(ReflectionList &reflections) {
    for (std::size_t i = 0; i < reflections.size(); ++i) {
      Reflection &r = reflections[i];
      int size_x = r.get_shoebox()[1] - r.get_shoebox()[0];
      int size_y = r.get_shoebox()[3] - r.get_shoebox()[2];
      int size_z = r.get_shoebox()[5] - r.get_shoebox()[4];
      r.set_image(flex_int(flex_grid<>(size_z, size_y, size_x)));
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
   */
  inline
  void copy_single_image_pixels(const flex_int &image, int array_index,
      const flex_int &index, ReflectionList &reflections) {
    for (std::size_t i = 0; i < index.size(); ++i) {

      // Get a reference to a reflection
      Reflection &r = reflections[index[i]];
      int6 shoebox = r.get_shoebox();
      int i0 = shoebox[0], i1 = shoebox[1];
      int j0 = shoebox[2], j1 = shoebox[3];
      int k0 = shoebox[4];
      int k = array_index - k0;

      // Get the image size
      flex_int::index_type image_size = index.accessor().all();

      // Readjust the area to loop over to ensure we're within image bounds
      int jj0 = j0 >= 0 ? j0 : 0;
      int ii0 = i0 >= 0 ? i0 : 0;
      int jj1 = j1 <= image_size[0] ? j1 : image_size[0];
      int ii1 = i1 <= image_size[1] ? i1 : image_size[1];

      // Get the reflection profile
      flex_int profile = r.get_image();

      // Copy the image pixels
      for (int jj = jj0; jj < jj1; ++jj) {
        for (int ii = ii0; ii < ii1; ++ii) {
          int j = jj - j0;
          int i = ii - i0;
          profile(k, j, i) = image(jj, ii);
        }
      }

      // Set the reflection profile
      r.set_image(profile);
    }
  }

}} // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_INTEGRATION_REFLECTION_PROFILE_HELPERS_H */
