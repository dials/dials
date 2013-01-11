
#ifndef DIALS_GEOMETRY_TRANSFORM_FROM_BEAM_VECTOR_TO_DETECTOR_H
#define DIALS_GEOMETRY_TRANSFORM_FROM_BEAM_VECTOR_TO_DETECTOR_H

#include <exception>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include "../detector_coordinate_system.h"

namespace dials { namespace geometry { namespace transform {

/** 
 * Class to represent a geometry transform from beam vector to detector
 * coordinates 
 */
class FromBeamVectorToDetector {

public:

    /** Default constructor */
    FromBeamVectorToDetector() : _distance(0.0) {}

    /** 
     * Initialise the transform from the detector coordinate system. The
     * detector coordinate system needs to be scaled in pixel units.
     * @param dcs The detector coordinate system
     * @param origin The origin of the detector coordinate system
     * @param distance The distance from the detector to the crystal
     */
    FromBeamVectorToDetector(DetectorCoordinateSystem dcs,
                             scitbx::vec2 <double> origin,
                             double distance) 
        : _x_axis(dcs.get_x_axis()),
          _y_axis(dcs.get_y_axis()),
          _normal(dcs.get_normal()),
          _origin(origin),
          _distance(distance) {}

public:

    /**
     * Apply the transform to a beam vector
     * @param s1 The beam vector
     * @returns A 2 element vector containing the pixel coordinates
     * @throws std::runtime_error if beam vector does not intersect with the
     *         detector plane
     */
    scitbx::vec2 <double> apply(scitbx::vec3 <double> s1) {
        double s1_dot_n = s1 * _normal;
        if (_distance * s1_dot_n <= 0) {
            throw std::runtime_error("Can't map beam vector to detector");
        }
        return scitbx::vec2 <double> (
            _origin[0] + _distance * (s1 * _x_axis) / s1_dot_n,
            _origin[1] + _distance * (s1 * _y_axis) / s1_dot_n);
    }

private:

    scitbx::vec3 <double> _x_axis;
    scitbx::vec3 <double> _y_axis;
    scitbx::vec3 <double> _normal;
    scitbx::vec2 <double> _origin;
    double _distance;
};

}}} // namespace = dials::geometry::transform

#endif // DIALS_GEOMETRY_TRANSFORM_FROM_BEAM_VECTOR_TO_DETECTOR_H
