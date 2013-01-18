
#ifndef DIALS_GEOMETRY_TRANSFORM_FROM_BEAM_VECTOR_TO_IMAGE_VOLUME_H
#define DIALS_GEOMETRY_TRANSFORM_FROM_BEAM_VECTOR_TO_IMAGE_VOLUME_H

#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include "../../equipment/goniometer.h"
#include "../../equipment/detector.h"
#include "from_beam_vector_to_detector.h"

namespace dials { namespace geometry { namespace transform {

/** A transform from beam vector to image volume */
class FromBeamVectorToImageVolume {

public:

    /** Default constructor */
    FromBeamVectorToImageVolume() {}

    /**
     * Initialise the transform from the detector and goniometer.
     * @param detector The detector struct
     * @param goniometer The goniometer struct
     */
    FromBeamVectorToImageVolume(equipment::Detector detector,
                                equipment::Goniometer goniometer)
        : from_beam_vector_to_detector_(
            DetectorCoordinateSystem(
                detector.get_x_axis(),
                detector.get_y_axis(), 
                detector.get_normal()),
                detector.get_pixel_size(),
                detector.get_origin(),
                detector.get_distance()),
          goniometer_(goniometer) {}

    /**
     * Apply the transform to a beam vector and rotation angle.
     * @param s1 The beam vector
     * @param phi The rotation angle
     * @returns The (x, y, z) image volume coordinate
     */
    scitbx::vec3 <double> apply(scitbx::vec3 <double> s1, double phi) {
        scitbx::vec2 <double> xy;
        try {
            xy = from_beam_vector_to_detector_.apply(s1);
        } catch(error) {
            xy = scitbx::vec2 <double> (-1, -1);
        }
        return scitbx::vec3 <double> (
            xy[0], xy[1],
            (goniometer_.get_frame_from_angle(phi) - 
             goniometer_.get_starting_frame()));
    }

    /**
     * Apply the transform to an array of beam vectors and rotation angles.
     * @param s1 The array of beam vectors
     * @param phi The array of rotation angles
     * @returns The array of image volume coordinates.
     */
    flex_vec3_double apply(flex_vec3_double s1, scitbx::af::flex_double phi) {
        DIALS_ASSERT(s1.size() == phi.size());
        flex_vec3_double result(s1.size());
        for (int i = 0; i < s1.size(); ++i) {
            result[i] = apply(s1[i], phi[i]);
        }
        return result;
    }

private:

    FromBeamVectorToDetector from_beam_vector_to_detector_;
    equipment::Goniometer goniometer_;
};

}}} // namespace dials::geometry::transform

#endif // DIALS_GEOMETRY_TRANSFORM_FROM_BEAM_VECTOR_TO_IMAGE_VOLUME_H