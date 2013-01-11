
#ifndef DIALS_GEOMETRY_TRANSFORM_FROM_DETECTOR_TO_BEAM_VECTOR_H
#define DIALS_GEOMETRY_TRANSFORM_FROM_DETECTOR_TO_BEAM_VECTOR_H

#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include "../detector_coordinate_system.h"

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
     * @param origin The origin of the detector coordinate system
     * @param distance The distance from the detector to the crystal
     */
    FromDetectorToBeamVector(DetectorCoordinateSystem dcs,
                             scitbx::vec2 <double> origin,
                             double distance) 
        : _x_axis(dcs.get_x_axis()),
          _y_axis(dcs.get_y_axis()),
          _distance_scaled_normal(dcs.get_normal() * distance),
          _origin(origin) {}

public:

    /**
     * Apply the transform to a beam vector
     * @param xy The detector coordinates
     * @returns The beam vector
     */
    scitbx::vec3 <double> apply(scitbx::vec2 <double> xy) {
        return (xy[0] - _origin[0]) * _x_axis + 
               (xy[1] - _origin[1]) * _y_axis + 
               _distance_scaled_normal;
    }

private:

    scitbx::vec3 <double> _x_axis;
    scitbx::vec3 <double> _y_axis;
    scitbx::vec3 <double> _distance_scaled_normal;
    scitbx::vec2 <double> _origin;
};

}}} // namespace = dials::geometry::transform

#endif // DIALS_GEOMETRY_TRANSFORM_FROM_DETECTOR_TO_BEAM_VECTOR_H
