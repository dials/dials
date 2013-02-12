
#ifndef DIALS_SPOT_PREDICTION_SUBTRACT_BACKGROUND_H
#define DIALS_SPOT_PREDICTION_SUBTRACT_BACKGROUND_H

#include <algorithm>
#include <scitbx/vec3.h>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/ref_reductions.h>
#include <scitbx/math/mean_and_variance.h>
#include "../error.h"
#include "../array_family/array_types.h"
#include "background_intensity.h"
#include "../reflection/reflection.h"

namespace dials { namespace integration {

/** A class to subtract the background intensity from the reflection profile */
class SubtractBackground {

public:

    /** Default constructor */
    SubtractBackground() {}

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
    SubtractBackground(const scitbx::af::flex_int &image_volume,
                       const scitbx::af::flex_int &reflection_mask,
                       double delta = 0.1,
                       double max_iter_frac = 0.1,
                       int min_pixels = 10)
        : image_volume_(image_volume),
          reflection_mask_(reflection_mask),
          background_intensity_(delta, max_iter_frac),
          min_pixels_(min_pixels)
    {
        DIALS_ASSERT(are_image_sizes_valid());
    }

    /**
     * Calculate the background intensity for a single reflection and subtract
     * it from the image pixel values.
     *
     * @todo In the XDS paper, the background intensity value is over estimated
     *       for strong reflections and in adjusted using the modelled
     *       intensity profile in the xds frame. This needs to be done.
     *
     * @param roi The region of interest
     */
    double subtract(int index, scitbx::af::tiny <int, 6> roi)
    {
        // Check the roi is valid
        DIALS_ASSERT(is_roi_valid(roi));

        // Number of pixels in the ROI
        int num_roi = (roi[1] - roi[0]) *
                      (roi[3] - roi[2]) *
                      (roi[5] - roi[4]);

        // Allocate memory for a temp array
        scitbx::af::flex_double data(num_roi);

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
        double background_value = background_intensity_.calculate(
                                    scitbx::af::flex_double_ref(
                                        data.begin(),
                                        data_index));

        // Loop through elements, subtract background and ensure >= 0
        for (int k = roi[4]; k < roi[5]; ++k) {
            for (int j = roi[2]; j < roi[3]; ++j) {
                for (int i = roi[0]; i < roi[1]; ++i) {
                    if (reflection_mask_(k, j, i) == index) {
                        image_volume_(k, j, i) -= background_value;
                        if (image_volume_(k, j, i) < 0) {
                            image_volume_(k, j, i) = 0;
                        }
                    }
                }
            }
        }

        return background_value;
    }

    /**
     * Subtract the background for all reflections
     * @param reflections The array of reflections
     * @returns The a boolean array containing the status for each reflection.
     *          True/False was the background successfully subtracted
     */
    scitbx::af::flex_bool subtract(ReflectionList &reflections) {
        scitbx::af::flex_bool result(reflections.size());
        for (int i = 0; i < reflections.size(); ++i) {
            try {
                reflections[i].set_background_intensity(
                    subtract(reflections[i].get_mask_index(),
                         reflections[i].get_region_of_interest()));
                result[i] = true;
            } catch(error) {
                result[i] = false;
            }
        }
        return result;
    }

    /**
     * If the pixel does not belong to a reflection then set the image volume
     * pixel value to the given value.
     * @param value The pixel value to set for non reflection pixels.
     */
    void set_non_reflection_value(int value) {
        for (int i = 0; i < image_volume_.size(); ++i) {
            if (reflection_mask_[i] == -1) {
                image_volume_[i] = value;
            }
        }
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
    bool is_roi_valid(scitbx::af::tiny <double, 6> roi) const {
        return roi[0] >= 0 && roi[1] <= image_volume_.accessor().all()[2] &&
               roi[2] >= 0 && roi[3] <= image_volume_.accessor().all()[1] &&
               roi[4] >= 0 && roi[5] <= image_volume_.accessor().all()[0];
    }

    scitbx::af::flex_int image_volume_;
    scitbx::af::flex_int reflection_mask_;
    scitbx::vec3 <int> roi_size_;
    BackgroundIntensity background_intensity_;
    int min_pixels_;
};

}} // namespace = dials::spot_prediction

#endif // DIALS_SPOT_PREDICTION_SUBTRACT_BACKGROUND_H
