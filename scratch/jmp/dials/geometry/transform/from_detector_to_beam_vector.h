
#ifndef DIALS_GEOMETRY_TRANSFORM_FROM_DETECTOR_TO_BEAM_VECTOR_H
#define DIALS_GEOMETRY_TRANSFORM_FROM_DETECTOR_TO_BEAM_VECTOR_H

#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include "../../equipment/detector.h"

namespace dials { namespace geometry { namespace transform {

/**
 * Class to represent a geometry transform from detector coordinates to beam
 * vector
 */
class FromDetectorToBeamVector {

public:

    /** Default constructor */
    FromDetectorToBeamVector() {}

    /**
     * Initialise the transform from the detector coordinate system. The
     * detector coordinate system needs to be scaled in mm.
     * @param dcs The detector coordinate system
     * @param pixel_size The size of the pixels in mm
     * @param origin The origin of the detector coordinate system
     * @param distance The distance from the detector to the crystal
     */
    FromDetectorToBeamVector(const equipment::Detector &detector)
        : x_axis_(detector.get_x_axis().normalize() *
                    detector.get_pixel_size()[0]),
          y_axis_(detector.get_y_axis().normalize() *
                    detector.get_pixel_size()[1]),
          distance_scaled_normal_(detector.get_normal().normalize() *
                    detector.get_distance()),
          origin_(detector.get_origin()) {}

public:

    /**
     * Apply the transform to a beam vector
     * @param xy The detector coordinates
     * @returns The beam vector
     */
    scitbx::vec3 <double> apply(scitbx::vec2 <double> xy) const {
        return (xy[0] - origin_[0]) * x_axis_ +
               (xy[1] - origin_[1]) * y_axis_ +
               distance_scaled_normal_;
    }

private:

    scitbx::vec3 <double> x_axis_;
    scitbx::vec3 <double> y_axis_;
    scitbx::vec3 <double> distance_scaled_normal_;
    scitbx::vec2 <double> origin_;
};

}}} // namespace = dials::geometry::transform

#endif // DIALS_GEOMETRY_TRANSFORM_FROM_DETECTOR_TO_BEAM_VECTOR_H
