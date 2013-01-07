
#ifndef DIALS_EQUIPMENT_DETECTOR_H
#define DIALS_EQUIPMENT_DETECTOR_H

#include <scitbx/vec2.h>
#include <scitbx/vec3.h>

namespace dials {
namespace equipment {

class detector {

public:

    detector() : _distance(0.0) {}
    
    detector(scitbx::vec3 <double> axis_x,
             scitbx::vec3 <double> axis_y,
             scitbx::vec3 <double> normal,
             scitbx::vec2 <double> origin,
             scitbx::vec2 <double> pixel_size,
             double distance)
        : _axis_x(axis_x),
          _axis_y(axis_y),
          _normal(normal),
          _origin(origin),
          _pixel_size(pixel_size),
          _distance(distance) {}

public:

    scitbx::vec3 <double> get_axis_x() {
        return _axis_x;
    }
    
    scitbx::vec3 <double> get_axis_y() {
        return _axis_y;
    }
    
    scitbx::vec3 <double> get_normal() {
        return _normal;
    }

    scitbx::vec2 <double> get_origin() {
        return _origin;
    }
    
    scitbx::vec2 <double> get_pixel_size() {
        return _pixel_size;
    }

    double get_distance() {
        return _distance;
    }
    
    void set_axis_x(scitbx::vec3 <double> axis_x) {
        _axis_x = axis_x;
    }
    
    void set_axis_y(scitbx::vec3 <double> axis_y) {
        _axis_y = axis_y;
    }
    
    void set_normal(scitbx::vec3 <double> normal) {
        _normal = normal;
    }
        
    void set_origin(scitbx::vec2 <double> origin) {
        _origin = origin;
    }
    
    void set_pixel_size(scitbx::vec2 <double> pixel_size) {
        _pixel_size = pixel_size;
    }
    
    void set_distance(double distance) {
        _distance = distance;
    }
        
private:

    scitbx::vec3 <double> _axis_x;
    scitbx::vec3 <double> _axis_y;
    scitbx::vec3 <double> _normal;
    scitbx::vec2 <double> _origin;
    scitbx::vec2 <double> _pixel_size;
    double _distance;
};

}
}

#endif // DIALS_EQUIPMENT_DETECTOR_H
