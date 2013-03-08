/*
 * subtract_background.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_INTEGRATION_SUBTRACT_BACKGROUND_H
#define DIALS_ALGORITHMS_INTEGRATION_SUBTRACT_BACKGROUND_H

#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/flex_types.h>
#include <dials/model/data/reflection.h>
#include <dials/error.h>
#include "background_intensity.h"

namespace dials { namespace algorithms {

  using scitbx::af::int6;
  using scitbx::af::flex_bool;
  using scitbx::af::flex_int;
  using scitbx::af::flex_double;
  using dials::model::Reflection;
  using dials::model::ReflectionList;

  /** A class to subtract the background intensity from the reflection profile */
  class SubtractBackground {

  public:

    /**
     * Initialise the class with parameters
     * @param image_volume The 3D image volume array
     * @param reflection_mask The 3D reflection mask
     * @param delta The deviation from normal
     * @param max_iter The maximum numner of iterations as a fraction of the
     *        elements in the input data array
     * @param min_pixels The minimum number of pixels needed to calculate
     *        the background intensity
     */
    SubtractBackground(const flex_int &image_volume,
                       const flex_int &reflection_mask,
                       int min_pixels = 10,
                       double n_sigma = -1)
      : image_volume_(image_volume),
        reflection_mask_(reflection_mask),
        min_pixels_(min_pixels),
        n_sigma_(n_sigma) {
      DIALS_ASSERT(are_image_sizes_valid());
    }

    /**
     * Calculate the background intensity for a single reflection and subtract
     * it from the image pixel values.
     *
     * @todo In the XDS paper, the background intensity value is over estimated
     *       for strong reflections and is adjusted using the modelled
     *       intensity profile in the xds frame. This needs to be done.
     *
     * @param roi The region of interest
     */
    double operator()(int index, int6 roi)
    {
      // Check the roi is valid
      DIALS_ASSERT(is_roi_valid(roi));

      // Number of pixels in the ROI
      int num_roi = (roi[1] - roi[0]) * (roi[3] - roi[2]) * (roi[5] - roi[4]);

      // Allocate memory for a temp array
      flex_double data(num_roi);

      // Copy the image pixels into a temp array
      int data_index = 0;
      for (int k = roi[4]; k < roi[5]; ++k) {
        for (int j = roi[2]; j < roi[3]; ++j) {
          for (int i = roi[0]; i < roi[1]; ++i) {
            if (reflection_mask_(k, j, i) == index) {
                data[data_index++] = image_volume_(k, j, i);
            }
          }
        }
      }

      // Ensure we have enough pixels to calculate the background
      DIALS_ASSERT(data_index > min_pixels_);

      // Calculate the background value
      double background_value = background_intensity(
        data, min_pixels_, n_sigma_);

      // Loop through elements, subtract background and ensure >= 0
      for (int k = roi[4]; k < roi[5]; ++k) {
        for (int j = roi[2]; j < roi[3]; ++j) {
          for (int i = roi[0]; i < roi[1]; ++i) {
            if (reflection_mask_(k, j, i) == index) {
              image_volume_(k, j, i) -= background_value;
            }
          }
        }
      }

      // Return the calculated background value
      return background_value;
    }

    /**
     * Subtract the background for all reflections
     * @param reflections The array of reflections
     * @returns The a boolean array containing the status for each reflection.
     *          True/False was the background successfully subtracted
     */
    flex_bool operator()(ReflectionList &reflections) {
      flex_bool result(reflections.size());
      for (int i = 0; i < reflections.size(); ++i) {
        try {
//          reflections[i].set_background_intensity(
//              subtract(reflections[i].get_mask_index(),
//                   reflections[i].get_shoebox()));
          result[i] = true;
        } catch(error) {
          result[i] = false;
        }
      }
      return result;
    }

  private:

    /** Ensure the images are 3D and of the same size */
    bool are_image_sizes_valid() const {
      return image_volume_.accessor().all().size() == 3
          && reflection_mask_.accessor().all().size() == 3
          && (image_volume_.accessor().all() ==
              reflection_mask_.accessor().all()).all_eq(true);
    }

    /** Check the roi is valid */
    bool is_roi_valid(int6 roi) const {
      return roi[0] >= 0 && roi[1] <= image_volume_.accessor().all()[2]
          && roi[2] >= 0 && roi[3] <= image_volume_.accessor().all()[1]
          && roi[4] >= 0 && roi[5] <= image_volume_.accessor().all()[0];
    }

    flex_int image_volume_;
    flex_int reflection_mask_;
    int min_pixels_;
    double n_sigma_;
  };

}} // namespace = dials::algorithms

#endif // DIALS_ALGORITHMS_INTEGRATION_SUBTRACT_BACKGROUND_H
