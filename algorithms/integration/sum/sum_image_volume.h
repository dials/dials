/*
 * sum_image_volume.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_INTEGRATION_SUM_SUM_IMAGE_VOLUME_H
#define DIALS_ALGORITHMS_INTEGRATION_SUM_SUM_IMAGE_VOLUME_H

#include <dials/array_family/reflection_table.h>
#include <dials/model/data/image_volume.h>
#include <dials/model/data/observation.h>

namespace dials { namespace algorithms {

  using dials::model::ImageVolume;
  using dials::model::Intensity;
  using dials::model::MultiPanelImageVolume;

  /**
   * Compute the summation intensity from a single reflection
   */
  template <typename FloatType>
  Intensity sum_image_volume(std::size_t index,
                             int6 bbox,
                             ImageVolume<FloatType> volume,
                             bool success) {
    // Trim the bbox
    int6 trimmed_bbox = volume.trim_bbox(bbox);

    // Do the summation
    Summation<FloatType> summation(
      volume.extract_data(trimmed_bbox).const_ref(),
      volume.extract_background(trimmed_bbox).const_ref(),
      volume.extract_mask(trimmed_bbox, index).const_ref());

    // Return the result
    Intensity result;
    result.observed.value = summation.intensity();
    result.observed.variance = summation.variance();
    result.background.value = summation.background();
    result.background.variance = summation.background_variance();
    result.observed.success = summation.success() && success;
    return result;
  }

  /**
   * Compute summation intensities from all reflections
   */
  template <typename FloatType>
  af::shared<Intensity> sum_multi_panel_image_volume(
    af::reflection_table reflections,
    MultiPanelImageVolume<FloatType> volume) {
    DIALS_ASSERT(reflections.contains("bbox"));
    DIALS_ASSERT(reflections.contains("panel"));
    DIALS_ASSERT(reflections.contains("fraction"));
    af::const_ref<int6> bbox = reflections["bbox"];
    af::const_ref<std::size_t> panel = reflections["panel"];
    af::const_ref<double> fraction = reflections["fraction"];
    af::shared<Intensity> intensity(bbox.size());
    for (std::size_t i = 0; i < bbox.size(); ++i) {
      intensity[i] =
        sum_image_volume(i, bbox[i], volume.get(panel[i]), fraction[i] > 1.0 - 1e-6);
    }
    return intensity;
  }

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_INTEGRATION_SUM_SUM_IMAGE_VOLUME_H
