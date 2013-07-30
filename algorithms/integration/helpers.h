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
#ifndef DIALS_ALGORITHMS_INTEGRATION_HELPERS_H
#define DIALS_ALGORITHMS_INTEGRATION_HELPERS_H

#include <limits>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/tiny_types.h>
#include <dials/model/data/reflection.h>
#include <dials/algorithms/image/threshold/unimodal.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using scitbx::af::small;
  using scitbx::af::int2;
  using scitbx::af::int6;
  using scitbx::af::flex_bool;
  using scitbx::af::flex_double;
  using dials::model::Reflection;
  using dials::model::ReflectionList;
  using dials::algorithms::maximum_deviation;

  /**
   * Check if the bounding box has points outside the image range.
   * @param bbox The bounding box
   * @param image_size The image size
   * @param scan_range The scan range
   * @returns True/False
   */
  bool is_bbox_outside_image_range(int6 bbox, small<long,10> image_size,
      int2 scan_range) {
    DIALS_ASSERT(image_size.size() == 2);
    return bbox[0] < 0 || bbox[1] >= image_size[1] ||
           bbox[2] < 0 || bbox[3] >= image_size[0] ||
           bbox[4] < scan_range[0] || bbox[5] >= scan_range[1];
  }

  /**
   * Check if the bounding box has points that cover bad pixels
   * @param bbox The bounding box
   * @param mask The mask array
   * @returns True/False
   */
  bool does_bbox_contain_bad_pixels(int6 bbox, const flex_bool &mask) {
    DIALS_ASSERT(mask.accessor().all().size() == 2);
    for (int j = bbox[2]; j < bbox[3]; ++j) {
      for (int i = bbox[0]; i < bbox[1]; ++i) {
        if (mask(j, i) == false) {
          return true;
        }
      }
    }
    return false;
  }

  /**
   * Check if the bounding box is valid in terms of the detector mask
   * @param bbox The bounding box
   * @param mask The mask array
   * @param scan_range The scan range
   * @returns True/False
   */
  bool is_bbox_valid(int6 bbox, const flex_bool &mask, int2 scan_range) {
    return !(is_bbox_outside_image_range(
        bbox, mask.accessor().all(), scan_range) ||
        does_bbox_contain_bad_pixels(bbox, mask));
  }

  /**
   * Filter the reflection based on the detector mask
   * @param reflection The reflection
   * @param mask The detector mask
   * @param scan_range The scan range
   */
  void filter_by_detector_mask(Reflection &reflection, const flex_bool &mask,
      int2 scan_range) {
    // get the bounding box
    int6 bbox = reflection.get_bounding_box();

    // Set whether the reflection is valid or not
    reflection.set_valid(is_bbox_valid(bbox, mask, scan_range));
  }

  /**
   * Filter the reflection list based on the detector mask
   * @param reflection The reflection
   * @param mask The detector mask
   * @param scan_range The scan range
   */
  void filter_by_detector_mask(ReflectionList &reflections,
      const flex_bool &mask, int2 scan_range) {
    for (std::size_t i = 0; i < reflections.size(); ++i) {
      filter_by_detector_mask(reflections[i], mask, scan_range);
    }
  }

}} // namespace dials::algorithms

#endif // DIALS_ALGORITHMS_INTEGRATION_HELPERS_H
