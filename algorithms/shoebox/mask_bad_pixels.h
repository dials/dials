/*
 * mask_bad_pixels.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_SHOEBOX_MASK_BAD_PIXELS_H
#define DIALS_ALGORITHMS_SHOEBOX_MASK_BAD_PIXELS_H

#include <scitbx/array_family/flex_types.h>
#include <dials/model/data/reflection.h>
#include <dials/algorithms/shoebox/mask_code.h>
#include <dials/error.h>

namespace dials { namespace algorithms { namespace shoebox {

  using scitbx::af::int6;
  using scitbx::af::flex_bool;
  using scitbx::af::flex_int;
  using dials::model::Reflection;
  using dials::model::ReflectionList;

  /**
   * A class to mask bad pixels from the reflection shoebox
   */
  class MaskBadPixels {
  public:

    /**
     * Initialise the algorithm with the detector mask
     * @param detector_mask The mask of good/bad detector pixels.
     */
    MaskBadPixels(const flex_bool &detector_mask)
      : detector_mask_(detector_mask),
        detector_size_(detector_mask.accessor().all()) {
      DIALS_ASSERT(detector_size_.size() == 2);
      DIALS_ASSERT(detector_size_[0] > 0 && detector_size_[1] > 0);
    }

    /**
     * Set all the shoebox mask values for the reflection to value
     * @param reflection The reflection
     */
    void operator()(Reflection &reflection) const {

      // Get the mask and roi from the reflection
      flex_int mask = reflection.get_shoebox_mask();
      flex_int::index_type size = mask.accessor().all();
      int6 roi = reflection.get_bounding_box();

      // Check sizes are ok
      DIALS_ASSERT(size.size() == 3);
      DIALS_ASSERT(size[2] == (roi[1] - roi[0]));
      DIALS_ASSERT(size[1] == (roi[3] - roi[2]));
      DIALS_ASSERT(size[0] == (roi[5] - roi[4]));

      // Loop through all the pixels in the reflection mask and check them
      // against the detector mask. If the detector mask pixel is ok then
      // set the reflection mask pixel to the given value, otherwise set it
      // to zero.
      for (std::size_t j = 0; j < size[1]; ++j) {
        for (std::size_t i = 0; i < size[2]; ++i) {
          int di = roi[0] + i;
          int dj = roi[2] + j;
          bool mask_value = false;
          if (di >= 0 && di < detector_size_[1] &&
              dj >= 0 && dj < detector_size_[0]) {
            mask_value = detector_mask_(dj, di);
          }
          for (std::size_t k = 0; k < size[0]; ++k) {
            mask(k, j, i) = (mask_value ? Valid : 0);
          }
        }
      }
    }

    /**
     * Set all the shoebox mask values for all the reflections to value
     * @param reflections The reflection list
     */
    void operator()(ReflectionList &reflections) const {
      for (std::size_t i = 0; i < reflections.size(); ++i) {
        if (reflections[i].is_valid()) {
          this->operator()(reflections[i]);
        }
      }
    }

  private:
    flex_bool detector_mask_;
    flex_bool::index_type detector_size_;
  };

}}} // namespace dials::algorithms::shoebox

#endif /* DIALS_ALGORITHMS_SHOEBOX_MASK_BAD_PIXELS_H */
