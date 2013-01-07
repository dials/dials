
#ifndef DIALS_GEOMETRY_TRANSFORM_FROM_BEAM_VECTOR_TO_DETECTOR_H
#define DIALS_GEOMETRY_TRANSFORM_FROM_BEAM_VECTOR_TO_DETECTOR_H

#include <exception>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>

namespace dials { namespace geometry { namespace transform {

class from_beam_vector_to_detector {

public:

    from_beam_vector_to_detector() : _distance(0.0) {}

    from_beam_vector_to_detector(scitbx::vec3 <double> axis_x,
                                 scitbx::vec3 <double> axis_y,
                                 scitbx::vec2 <double> origin,
                                 double distance)
        : _axis_x(axis_x),
          _axis_y(axis_y),
          _normal(axis_x.cross(axis_y).normalize()),
          _origin(origin),
          _distance(distance) {}

    from_beam_vector_to_detector(scitbx::vec3 <double> axis_x,
                                 scitbx::vec3 <double> axis_y,
                                 scitbx::vec3 <double> normal,
                                 scitbx::vec2 <double> origin,
                                 double distance) 
        : _axis_x(axis_x),
          _axis_y(axis_y),
          _normal(normal),
          _origin(origin),
          _distance(distance) {}

public:

    scitbx::vec2 <double> apply(scitbx::vec3 <double> s1) {

        double s1_dot_n = s1 * _normal;
        if (_distance * s1_dot_n <= 0) {
            throw std::runtime_error("Can't map beam vector to detector");
        }
        
        return scitbx::vec2 <double> (
            _origin[0] + _distance * (s1 * _axis_x) / s1_dot_n,
            _origin[1] + _distance * (s1 * _axis_y) / s1_dot_n);
    }

private:

    scitbx::vec3 <double> _axis_x;
    scitbx::vec3 <double> _axis_y;
    scitbx::vec3 <double> _normal;
    scitbx::vec2 <double> _origin;
    double _distance;
};

}}}

#endif // DIALS_GEOMETRY_TRANSFORM_FROM_BEAM_VECTOR_TO_DETECTOR_H
