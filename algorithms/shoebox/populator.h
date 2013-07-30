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
  using dials::model::Reflection;
  using dials::model::ReflectionList;

  class Populator {
  public:

    Populator(ReflectionList &reflections, const flex_bool &mask,
        const flex_double &gain_map, const flex_double &dark_map)
      : reflections_(reflections),
        mask_(mask),
        gain_map_(gain_map),
        dark_map_(dark_map) {
      allocate();
      initialize();
    }

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
          int k0 = bounding_box[4];
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

    void allocate() {
      // Allocate all the reflection profiles
      for (std::size_t i = 0; i < reflections_.size(); ++i) {
        Reflection &r = reflections_[i];
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

    void deallocate() {
      // Delete all the reflection profiles from memory
      for (std::size_t i = 0; i < reflections_.size(); ++i) {
        reflections_[i].set_shoebox(flex_double());
        reflections_[i].set_shoebox_mask(flex_int());
        reflections_[i].set_shoebox_background(flex_double());
        reflections_[i].set_transformed_shoebox(flex_double());
      }
    }

    void initialize() {
      store_frame_indices();
      mask_overlapping();
      mask_foreground();
    }

  private:

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

    void mask_overlapping() {

    }

    void mask_foreground() {

    }

    ReflectionList reflections_;
    flex_bool mask_;
    flex_double gain_map_;
    flex_double dark_map_;
    unordered_map<int, shared<int> > index_;
  };

}}} // namespace dials::algorithms::shoebox

#endif /* DIALS_ALGORITHMS_SHOEBOX_POPULATOR_H */
