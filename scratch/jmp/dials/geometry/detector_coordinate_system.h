
#ifndef DIALS_GEOMETRY_DETECTOR_COORDINATE_SYSTEM_H
#define DIALS_GEOMETRY_DETECTOR_COORDINATE_SYSTEM_H

#include <scitbx/vec2.h>
#include <scitbx/vec3.h>

namespace dials { namespace geometry {

class detector_coordinate_system {

public:

    detector_coordinate_system() {}

    detector_coordinate_system(scitbx::vec3 <double> axis_x,
                               scitbx::vec3 <double> axis_y)
        : _axis_x(axis_x),
          _axis_y(axis_y),
          _normal(axis_x.cross(axis_y)) {}
    
    detector_coordinate_system(scitbx::vec3 <double> axis_x,
                               scitbx::vec3 <double> axis_y,
                               scitbx::vec3 <double> normal)
        : _axis_x(axis_x),
          _axis_y(axis_y),
          _normal(normal) {}

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
    
    void set_axis_x(scitbx::vec3 <double> axis_x) {
        _axis_x = axis_x;
    }
    
    void set_axis_y(scitbx::vec3 <double> axis_y) {
        _axis_y = axis_y;
    }
    
    void set_normal(scitbx::vec3 <double> normal) {
        _normal = normal;
    }

public:

    static detector_coordinate_system in_pixel_units(
            detector_coordinate_system cs, scitbx::vec2 <double> pixel_size) {
        return detector_coordinate_system(
            cs.get_axis_x().normalize() / pixel_size[0],
            cs.get_axis_y().normalize() / pixel_size[1],
            cs.get_normal());
    }
    
    static detector_coordinate_system in_si_units(
            detector_coordinate_system cs, scitbx::vec2 <double> pixel_size) {
        return detector_coordinate_system(
            cs.get_axis_x().normalize() * pixel_size[0],
            cs.get_axis_y().normalize() * pixel_size[1],
            cs.get_normal());
    }

private:

    scitbx::vec3 <double> _axis_x;
    scitbx::vec3 <double> _axis_y;
    scitbx::vec3 <double> _normal;
};

}}

#endif // DIALS_GEOMETRY_DETECTOR_COORDINATE_SYSTEM_H
