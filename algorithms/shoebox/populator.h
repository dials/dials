/*
 * populator.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_SHOEBOX_POPULATOR_H
#define DIALS_ALGORITHMS_SHOEBOX_POPULATOR_H

#include <boost/unordered_map.hpp>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/shared.h>
#include <dials/model/data/reflection.h>
#include <dials/algorithms/shoebox/mask_code.h>
#include <dials/error.h>

namespace dials { namespace algorithms { namespace shoebox {

  using boost::unordered_map;
  using scitbx::af::int2;
  using scitbx::af::int6;
  using scitbx::af::flex_bool;
  using scitbx::af::flex_double;
  using scitbx::af::flex_int;
  using scitbx::af::flex_grid;
  using scitbx::af::shared;
  using scitbx::af::const_ref;
  using dials::model::Reflection;
  using dials::model::ReflectionList;

  /**
   * Class to allocate and populate reflection profiles from image data
   */
  class Populator {
  public:

    /**
     * Initialise the profiles.
     * @param reflections The list of reflections
     * @param mask The detector mask
     * @param gain_map The gain map
     * @param dark_map The dark map
     */
    Populator(ReflectionList &reflections, const flex_bool &mask,
        const flex_double &gain_map, const flex_double &dark_map)
      : reflections_(reflections),
        mask_(mask),
        gain_map_(gain_map),
        dark_map_(dark_map) {
      store_frame_indices();
    }

    /**
     * Add an image to the reflection shoebox
     * @param image The image pixels to add
     * @param image_index The index of the image
     */
    void add_image(flex_int &image, std::size_t image_index) {

      // Get the indices for this frame
      shared<int> indices = index_[image_index];

      // Loop through all the indices for this frame
      for (std::size_t i = 0; i < indices.size(); ++i) {

        // Get a reference to a reflection
        Reflection &r = reflections_[indices[i]];
        if (r.is_valid()) {
          int6 bounding_box = r.get_bounding_box();
          int i0 = bounding_box[0], i1 = bounding_box[1];
          int j0 = bounding_box[2], j1 = bounding_box[3];
          int k0 = bounding_box[4], k1 = bounding_box[5];
          int k = image_index - k0;

          // Get the image size
          flex_int::index_type image_size = image.accessor().all();

          // Readjust the area to loop over to ensure we're within image bounds
          int jj0 = j0 >= 0 ? j0 : 0;
          int ii0 = i0 >= 0 ? i0 : 0;
          int jj1 = j1 <= image_size[0] ? j1 : image_size[0];
          int ii1 = i1 <= image_size[1] ? i1 : image_size[1];

          // Get the reflection profile
          flex_double profile = r.get_shoebox();
          DIALS_ASSERT(profile.accessor().all()[0] == (k1 - k0));
          DIALS_ASSERT(profile.accessor().all()[1] == (j1 - j0));
          DIALS_ASSERT(profile.accessor().all()[2] == (i1 - i0));

          // Copy the image pixels
          for (int jj = jj0; jj < jj1; ++jj) {
            for (int ii = ii0; ii < ii1; ++ii) {
              int j = jj - j0;
              int i = ii - i0;
              profile(k, j, i) = gain_map_(jj, ii) * (
                  image(jj, ii) - dark_map_(jj, ii));
            }
          }

          // Set the reflection profile
          r.set_shoebox(profile);
        }
      }
    }

    /**
     * Get a mask composed of the shoeboxes on the image
     * @param image_index The index of the image
     * @param kernel_size The size to expand around the mask
     * @returns A mask for the image
     */
    flex_bool image_mask(int image_index, int2 kernel_size) {

      // Create the resulting mask
      flex_bool result(mask_.accessor(), false);
      flex_bool::index_type mask_size = result.accessor().all();
      shared<int> indices = index_[image_index];


      // Set all the shoebox pixels
      for (std::size_t i = 0; i < indices.size(); ++i) {

        // Get a reference to a reflection
        const Reflection &r = reflections_[indices[i]];
        if (r.is_valid()) {
          int6 bounding_box = r.get_bounding_box();
          int i0 = bounding_box[0] - kernel_size[0];
          int i1 = bounding_box[1] + kernel_size[0];
          int j0 = bounding_box[2] - kernel_size[1];
          int j1 = bounding_box[3] + kernel_size[1];

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
        result[i] = result[i] && mask_[i];
      }

      // Return the resulting mask
      return result;
    }

    /**
     * Get the reflection indices at the given image index
     * @param image_index The image index
     * @returns The reflection indices recorded on the image.
     */
    const shared<int> indices(int image_index) const {
      const const_ref<int> indices = index_.at(image_index).const_ref();
      shared<int> result(indices.size());
      for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] = indices[i];
      }
      return result;
    }

  private:

    /**
     * Get the image indices at which each reflection is recorded.
     */
    void store_frame_indices() {
      // For each reflection, Find the frames which it spans and copy an
      // index into the frame -> reflection list
      for (std::size_t i = 0; i < reflections_.size(); ++i) {
        Reflection& r = reflections_[i];
        if (r.is_valid()) {
          int f0 = r.get_bounding_box()[4];
          int f1 = r.get_bounding_box()[5];
          for (int f = f0; f < f1; ++f) {
            index_[f].push_back(i);
          }
        }
      }
    }

    ReflectionList reflections_;
    flex_bool mask_;
    flex_double gain_map_;
    flex_double dark_map_;
    unordered_map<int, shared<int> > index_;
  };

}}} // namespace dials::algorithms::shoebox

#endif /* DIALS_ALGORITHMS_SHOEBOX_POPULATOR_H */
