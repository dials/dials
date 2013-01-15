
#ifndef DIALS_SPOT_PREDICTION_SUBTRACT_BACKGROUND_H
#define DIALS_SPOT_PREDICTION_SUBTRACT_BACKGROUND_H

#include <algorithm>
#include <scitbx/vec3.h>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/math/mean_and_variance.h>
#include "background_intensity.h"

namespace dials { namespace spot_prediction {

/** flex array of vec3 doubles */
typedef scitbx::af::flex <scitbx::vec3 <double> >::type flex_vec3_double;

/** A class to subtract the background intensity from the reflection profile */
class SubtractBackground {

public:

    /** Default constructor */
    SubtractBackground() {}

    /**
     * Initialise the class with parameters
     * @param delta The deviation from normal
     * @param max_iter The maximum numner of iterations as a fraction of the
     *                 elements in the input data array
     */
    SubtractBackground(scitbx::af::flex_int image_volume,
                       scitbx::vec3 <int> roi_size,
                       double delta = 0.1, 
                       double max_iter_frac = 0.1) 
        : image_volume_(image_volume),
          roi_size_(roi_size),
          background_intensity_(delta, max_iter_frac) 
    {
        if (image_volume_.accessor().all().size() != 3) {
            throw std::runtime_error("Image volume must be 3D");
        }
    }

    /**
     * Calculate the background intensity
     * @param data The pixel array data
     */
    void subtract(flex_vec3_double xyz) 
    {
        // Allocate memory for a temp array
        scitbx::af::flex_double data((2 * roi_size_[0] + 1) * 
                                     (2 * roi_size_[1] + 1) *
                                     (2 * roi_size_[2] + 1));
        
        // Calculate the image stride
        int stride_x = image_volume_.accessor().all()[2];
        int stride_y = image_volume_.accessor().all()[1] * stride_x;

        // Loop through all the given reflection detector points.
        for (int index = 0; index < xyz.size(); ++index) {
            
            // Copy the image pixels into a temp array
            int x0 = (int)xyz[index][0] - roi_size_[2];
            int x1 = (int)xyz[index][0] + roi_size_[2];
            int y0 = (int)xyz[index][1] - roi_size_[1];
            int y1 = (int)xyz[index][1] + roi_size_[1];
            int z0 = (int)xyz[index][2] - roi_size_[0];
            int z1 = (int)xyz[index][2] + roi_size_[0];
            int data_index = 0;
            for (int k = z0; k <= z1; ++k) {
                for (int j = y0; j <= y1; ++j) {
                    for (int i = x0; i <= x1; ++i) {
                        int image_index = i + j * stride_x + k * stride_y;
                        data[data_index++] = image_volume_[image_index];
                    }
                }
            }
            
            // Calculate the background value
            double background = background_intensity_.calculate(
                                    scitbx::af::flex_double_const_ref(
                                        data.begin(),
                                        data_index));
        
            // Loop through elements, subtract background and ensure >= 0
            for (int k = z0; k <= z1; ++k) {
                for (int j = y0; j <= y1; ++j) {
                    for (int i = x0; i <= x1; ++i) {
                        int image_index = i + j * stride_x + k * stride_y;
                        image_volume_[image_index] -= background;
                        if (image_volume_[image_index] < 0) {
                            image_volume_[image_index] = 0;
                        }
                    }
                }
            }
        }
    }



private:

    scitbx::af::flex_int image_volume_;
    scitbx::vec3 <int> roi_size_;
    BackgroundIntensity background_intensity_;
};

}} // namespace = dials::spot_prediction

#endif // DIALS_SPOT_PREDICTION_SUBTRACT_BACKGROUND_H
