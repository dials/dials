

#ifndef DIALS_SPOT_PREDICTION_ROTATION_ANGLES_H
#define DIALS_SPOT_PREDICTION_ROTATION_ANGLES_H

#include <scitbx/constants.h>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <scitbx/array_family/flex_types.h>
#include <cctbx/miller.h>
#include <rstbx/diffraction/ewald_sphere.h>
#include "angle_filter.h"

namespace dials { namespace spot_prediction {

typedef cctbx::miller::index <> miller_index;
typedef scitbx::af::flex <miller_index>::type flex_miller_index;
typedef scitbx::af::flex_const_ref <miller_index>::type flex_miller_index_const_ref;

/**
 * A helper class that takes an array of miller indices and calculates the
 * set of valid rotation angles within a given range.
 */
class RotationAngles {

public:

    /**
     * Initialise the class and setup the underlying rotation angle class that
     * will calculate our rotation angles.
     * @param d_min The resolution
     * @param ub_matrix The UB matrix 
     * @param wavelength The incident beam wavelength
     * @param rotation_axis The crystal rotation axis
     * @param deg angles in degrees True/False
     */
    RotationAngles(double d_min, 
                   scitbx::mat3 <double> ub_matrix,
                   double wavelength, 
                   scitbx::vec3 <double> rotation_axis)
        : rotation_angles_(d_min, ub_matrix, wavelength, rotation_axis){}

    /**
     * Calculate the valid rotation angles from the given list of miller
     * indices. I.e Calculate the intersection angles and then check they are
     * within the given range of allowed rotation angles.
     * @param miller_indices The array of miller indices
     * @returns The array of rotation angles.
     */
    scitbx::af::flex_double calculate(flex_miller_index miller_indices,
                                          scitbx::af::flex_bool &status)
    {
//        miller_indices_ = flex_miller_index_const_ref(miller_indices.begin(), miller_indices.size());
//        array_indices_.resize(0);
        scitbx::af::flex_double result(scitbx::af::flex_grid <> (2, miller_indices.size()));
        status.resize(miller_indices.size());
        for (int i = 0; i < miller_indices.size(); ++i) {
            if (rotation_angles_(to_vec3_double(miller_indices[i]))) {
				scitbx::vec2 <double> angles = mod_2pi(
					rotation_angles_.get_intersection_angles());
                result(0, i) = angles[0];
                result(1, i) = angles[1];
                status[i] = true;
//                result.push_back(angles[0]);
//                array_indices_.push_back(i);
//                result.push_back(angles[1]);
//                array_indices_.push_back(i);
            } else {
                status[i] = false;
            }
        }
        return result;
    }
    
    /** Get the list of array indices corresponding the rotation angles */
    scitbx::af::shared <int> get_array_indices() {
        return array_indices_;
    }
    
    /** Get the miller indices corresponding to the rotation angles */
    flex_miller_index get_miller_indices() {
        flex_miller_index result(array_indices_.size());
        for (int i = 0; i < array_indices_.size(); ++i) {
            result[i] = miller_indices_[array_indices_[i]];
        }
        return result;
    }

private:

    /** Convert a miller::index struct to vec3 <double> */
    scitbx::vec3 <double> to_vec3_double(cctbx::miller::index <> i) {
        return scitbx::vec3 <double> ((double)i[0], (double)i[1], (double)i[2]);
    }

    rstbx::rotation_angles rotation_angles_;
    scitbx::af::shared <int> array_indices_;
    flex_miller_index_const_ref miller_indices_;
};

}} // namespace dials::spot_prediction

#endif // DIALS_SPOT_PREDICTION_ROTATION_ANGLES_H
