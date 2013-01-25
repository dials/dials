
#ifndef DIALS_EQUIPMENT_GONIOMETER_H
#define DIALS_EQUIPMENT_GONIOMETER_H

#include <scitbx/vec3.h>

namespace dials { namespace equipment {

/** The goniometer */
class Goniometer {

public:

    /** The default constructor */
    Goniometer() 
        : starting_angle_(0.0),
          oscillation_range_(0.0), 
          starting_frame_(0),
          num_frames_(0) {}


    /** 
     * Initialise the goniometer 
     * @param rotation_axis The goniometer rotation axis
     * @param starting_angle The starting rotation angle
     * @param oscillation_range The angular difference between successive frames
     * @param starting_frame The starting frame number
     * @param num_frames The number of frames
     */
    Goniometer(scitbx::vec3 <double> rotation_axis, 
               double starting_angle, 
               double oscillation_range,
               int starting_frame,
               int num_frames)
        : rotation_axis_(rotation_axis),
          starting_angle_(starting_angle),
          oscillation_range_(oscillation_range),
          starting_frame_(starting_frame),
          num_frames_(num_frames) {}
          
public:

    /** Get the rotation axis */
    scitbx::vec3 <double> get_rotation_axis() const {
        return rotation_axis_;
    }
   
    /** Get the starting angle */
    double get_starting_angle() const {
        return starting_angle_;
    }
    
    /** Get the oscillation range */
    double get_oscillation_range() const {
        return oscillation_range_;
    }

    /** Get the starting frame */
    int get_starting_frame() const {
        return starting_frame_;
    }
    
    /** Get the number of frames */
    int get_num_frames() const {
        return num_frames_;
    }
    
    /** Set the rotation axis */
    void set_rotation_axis(scitbx::vec3 <double> rotation_axis) {
        rotation_axis_ = rotation_axis;
    }
   
    /** Set the starting angle */
    void set_starting_angle(double starting_angle) {
        starting_angle_ = starting_angle;
    }
    
    /** Set the oscillation range */
    void set_oscillation_range(double oscillation_range) {
        oscillation_range_ = oscillation_range;
    }

    /** Set the starting frame */
    void set_starting_frame(int starting_frame) {
        starting_frame_ = starting_frame;
    }
    
    /** Set the number of frames */
    void set_num_frames(int num_frames) {
        num_frames_ = num_frames;
    }

public:

    /** 
     * Get the angle from the given frame 
     * @param frame The frame number
     * @returns The angle corresponding to the frame number
     */
    double get_angle_from_frame(double frame) {
        return starting_angle_ + (frame - starting_frame_) * oscillation_range_;
    }
    
    /**
     * Get the frame from the given angle
     * @param angle The angle of rotation
     * @returns The frame number corresponding to the rotation angle
     */
    double get_frame_from_angle(double angle) {
        return starting_frame_ + (angle - starting_angle_) / oscillation_range_;
    }
    
private:

    scitbx::vec3 <double> rotation_axis_;
    double starting_angle_;
    double oscillation_range_;
    int starting_frame_;
    int num_frames_;
};

}} // namespace = dials::equipment

#endif // DIALS_EQUIPMENT_GONIOMETER_H
