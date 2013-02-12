
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
    FromBeamVectorToImageVolume(const equipment::Detector &detector,
                                const equipment::Goniometer &goniometer)
        : from_beam_vector_to_detector_(detector),
          goniometer_(goniometer) {}

    /**
     * Apply the transform to a beam vector and rotation angle.
     * @param s1 The beam vector
     * @param phi The rotation angle
     * @returns The (x, y, z) image volume coordinate
     */
    scitbx::vec3 <double> apply(scitbx::vec3 <double> s1, double phi) const {
        scitbx::vec2 <double> xy = from_beam_vector_to_detector_.apply(s1);
        return scitbx::vec3 <double> (
            xy[0], xy[1],
            (goniometer_.get_zero_based_frame_from_angle(phi, true)));
    }

    /**
     * Apply the transform to an array of beam vectors and rotation angles.
     * @param s1 The array of beam vectors
     * @param phi The array of rotation angles
     * @param status The status array
     * @returns The array of image volume coordinates.
     */
    flex_vec3_double apply(const flex_vec3_double &s1,
                           const scitbx::af::flex_double &phi,
                           scitbx::af::flex_bool &status) const {
        DIALS_ASSERT(s1.size() == phi.size());
        status.resize(s1.size());
        flex_vec3_double result(s1.size());
        for (int i = 0; i < s1.size(); ++i) {
            try {
                result[i] = apply(s1[i], phi[i]);
                status[i] = true;
            } catch (error) {
                result[i] = scitbx::vec3 <double> (0, 0, 0);
                status[i] = false;
            }
        }
        return result;
    }

private:

    FromBeamVectorToDetector from_beam_vector_to_detector_;
    equipment::Goniometer goniometer_;
};

}}} // namespace dials::geometry::transform

#endif // DIALS_GEOMETRY_TRANSFORM_FROM_BEAM_VECTOR_TO_IMAGE_VOLUME_H
