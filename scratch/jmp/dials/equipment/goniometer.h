


#ifndef DIALS_EQUIPMENT_GONIOMETER_H
#define DIALS_EQUIPMENT_GONIOMETER_H

#include <scitbx/vec3.h>

namespace dials {
namespace equipment {

class goniometer {

public:

    goniometer() 
        : _starting_frame(0),
          _starting_angle(0.0),
          _oscillation_range(0.0) {}


    goniometer(scitbx::vec3 <double> rotation_axis, int starting_frame, 
               double starting_angle, double oscillation_range)
        : _rotation_axis(rotation_axis),
          _starting_frame(starting_frame),
          _starting_angle(starting_angle),
          _oscillation_range(oscillation_range) {}

public:

    scitbx::vec3 <double> get_rotation_axis() {
        return _rotation_axis;
    }

    int get_starting_frame() {
        return _starting_frame;
    }
    
    double get_starting_angle() {
        return _starting_angle;
    }
    
    double get_oscillation_range() {
        return _oscillation_range;
    }
    
    void set_rotation_axis(scitbx::vec3 <double> rotation_axis) {
        _rotation_axis = rotation_axis;
    }
    
    void set_starting_frame(int starting_frame) {
        _starting_frame = starting_frame;
    }
    
    void set_starting_angle(double starting_angle) {
        _starting_angle = starting_angle;
    }
    
    void set_oscillation_range(double oscillation_range) {
        _oscillation_range = oscillation_range;
    }

public:

    double get_angle_from_frame(double frame) {
        return _starting_angle + (frame - _starting_frame) * _oscillation_range;
    }
    
    double get_frame_from_angle(double angle) {
        return _starting_frame + (angle - _starting_angle) / _oscillation_range;
    }

private:

    scitbx::vec3 <double> _rotation_axis;
    int _starting_frame;
    double _starting_angle;
    double _oscillation_range;
};

}
}

#endif // DIALS_EQUIPMENT_GONIOMETER_H
