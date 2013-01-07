
#ifndef DIALS_EQUIPMENT_GONIOMETER_H
#define DIALS_EQUIPMENT_GONIOMETER_H

#include <scitbx/vec3.h>

namespace dials { namespace equipment {

/** The goniometer */
class goniometer {

public:

    /** The default constructor */
    goniometer() 
        : _starting_angle(0.0),
          _oscillation_range(0.0), 
          _starting_frame(0) {}


    /** 
     * Initialise the goniometer 
     * @param rotation_axis The goniometer rotation axis
     * @param starting_angle The starting rotation angle
     * @param oscillation_range The angular difference between successive frames
     * @param starting_frame The starting frame number
     */
    goniometer(scitbx::vec3 <double> rotation_axis, 
               double starting_angle, 
               double oscillation_range,
               int starting_frame)
        : _rotation_axis(rotation_axis),
          _starting_angle(starting_angle),
          _oscillation_range(oscillation_range),
          _starting_frame(starting_frame) {}
          
public:

    /** Get the rotation axis */
    scitbx::vec3 <double> get_rotation_axis() {
        return _rotation_axis;
    }
   
    /** Get the starting angle */
    double get_starting_angle() {
        return _starting_angle;
    }
    
    /** Get the oscillation range */
    double get_oscillation_range() {
        return _oscillation_range;
    }

    /** Get the starting frame */
    int get_starting_frame() {
        return _starting_frame;
    }
    
    /** Set the rotation axis */
    void set_rotation_axis(scitbx::vec3 <double> rotation_axis) {
        _rotation_axis = rotation_axis;
    }
   
    /** Set the starting angle */
    void set_starting_angle(double starting_angle) {
        _starting_angle = starting_angle;
    }
    
    /** Set the oscillation range */
    void set_oscillation_range(double oscillation_range) {
        _oscillation_range = oscillation_range;
    }

   /** Set the starting frame */
    void set_starting_frame(int starting_frame) {
        _starting_frame = starting_frame;
    }

public:

    /** 
     * Get the angle from the given frame 
     * @param frame The frame number
     * @returns The angle corresponding to the frame number
     */
    double get_angle_from_frame(double frame) {
        return _starting_angle + (frame - _starting_frame) * _oscillation_range;
    }
    
    /**
     * Get the frame from the given angle
     * @param angle The angle of rotation
     * @returns The frame number corresponding to the rotation angle
     */
    double get_frame_from_angle(double angle) {
        return _starting_frame + (angle - _starting_angle) / _oscillation_range;
    }

private:

    scitbx::vec3 <double> _rotation_axis;
    double _starting_angle;
    double _oscillation_range;
    int _starting_frame;
};

}
}

#endif // DIALS_EQUIPMENT_GONIOMETER_H
