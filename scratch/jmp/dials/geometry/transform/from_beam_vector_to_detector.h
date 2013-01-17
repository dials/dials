
#ifndef DIALS_GEOMETRY_TRANSFORM_FROM_BEAM_VECTOR_TO_DETECTOR_H
#define DIALS_GEOMETRY_TRANSFORM_FROM_BEAM_VECTOR_TO_DETECTOR_H

#include <exception>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include "../../error.h"
#include "../detector_coordinate_system.h"

namespace dials { namespace geometry { namespace transform {

/** 
 * Class to represent a geometry transform from beam vector to detector
 * coordinates 
 */
class FromBeamVectorToDetector {

public:

    /** Default constructor */
    FromBeamVectorToDetector() : distance_(0.0) {}

    /** 
     * Initialise the transform from the detector coordinate system. The
     * detector coordinate system needs to be scaled in pixel units.
     * @param dcs The detector coordinate system
     * @param pixel_size The size of the pixels in mm
     * @param origin The origin of the detector coordinate system
     * @param distance The distance from the detector to the crystal
     */
    FromBeamVectorToDetector(DetectorCoordinateSystem dcs,
                             scitbx::vec2 <double> pixel_size,
                             scitbx::vec2 <double> origin,
                             double distance) 
        : x_axis_(dcs.get_x_axis().normalize() / pixel_size[0]),
          y_axis_(dcs.get_y_axis().normalize() / pixel_size[1]),
          normal_(dcs.get_normal().normalize()),
          origin_(origin),
          distance_(distance) {}

public:

    /**
     * Apply the transform to a beam vector
     * @param s1 The beam vector
     * @returns A 2 element vector containing the pixel coordinates
     * @throws std::runtime_error if beam vector does not intersect with the
     *         detector plane
     */
    scitbx::vec2 <double> apply(scitbx::vec3 <double> s1) {
        double s1_dot_n = s1 * normal_;
        DIALS_ASSERT(distance_ * s1_dot_n > 0);
        return scitbx::vec2 <double> (
            origin_[0] + distance_ * (s1 * x_axis_) / s1_dot_n,
            origin_[1] + distance_ * (s1 * y_axis_) / s1_dot_n);
    }

private:

    scitbx::vec3 <double> x_axis_;
    scitbx::vec3 <double> y_axis_;
    scitbx::vec3 <double> normal_;
    scitbx::vec2 <double> origin_;
    double distance_;
};

}}} // namespace = dials::geometry::transform

#endif // DIALS_GEOMETRY_TRANSFORM_FROM_BEAM_VECTOR_TO_DETECTOR_H
