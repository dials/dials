
#ifndef DIALS_EQUIPMENT_DETECTOR_H
#define DIALS_EQUIPMENT_DETECTOR_H

#include <scitbx/vec2.h>
#include <scitbx/vec3.h>

namespace dials { namespace equipment {

/** The detector */
class Detector {

public:

    /** The default constructor */
    Detector() : distance_(0.0) {}
    
    /** 
     * Initialise the detector.
     * @param x_axis The x axis of the detector
     * @param y_axis The y axis of the detector
     * @param normal The detector normal
     * @param origin The detector origin (in pixels)
     * @param pixel_size The size of the individual pixels
     * @param size The size of the detector (in pixels)
     * @param distance The distance from the detector to the crystal origin
     */
    Detector(scitbx::vec3 <double> x_axis,
             scitbx::vec3 <double> y_axis,
             scitbx::vec3 <double> normal,
             scitbx::vec2 <double> origin,
             scitbx::vec2 <double> pixel_size,
             scitbx::vec2 <int> size,
             double distance)
        : x_axis_(x_axis),
          y_axis_(y_axis),
          normal_(normal),
          origin_(origin),
          pixel_size_(pixel_size),
          size_(size),
          distance_(distance) {}

public:

    /** Get the x axis */
    scitbx::vec3 <double> get_x_axis() const {
        return x_axis_;
    }
    
    /** Get the y axis */
    scitbx::vec3 <double> get_y_axis() const {
        return y_axis_;
    }
    
    /** Get the normal */
    scitbx::vec3 <double> get_normal() const {
        return normal_;
    }

    /** Get the detector origin */
    scitbx::vec2 <double> get_origin() const {
        return origin_;
    }
    
    /** Get the pixel size */
    scitbx::vec2 <double> get_pixel_size() const {
        return pixel_size_;
    }

    /** Get the detector size */
    scitbx::vec2 <int> get_size() const {
        return size_;
    }
    
    /** Get the distance from the crystal */
    double get_distance() const {
        return distance_;
    }
    
    /** Set the x axis */
    void set_x_axis(scitbx::vec3 <double> x_axis) {
        x_axis_ = x_axis;
    }
    
    /** Set the y axis */
    void set_y_axis(scitbx::vec3 <double> y_axis) {
        y_axis_ = y_axis;
    }
    
    /** Set the normal */
    void set_normal(scitbx::vec3 <double> normal) {
        normal_ = normal;
    }
        
    /** Set the origin */
    void set_origin(scitbx::vec2 <double> origin) {
        origin_ = origin;
    }
    
    /** Set the pixel size */
    void set_pixel_size(scitbx::vec2 <double> pixel_size) {
        pixel_size_ = pixel_size;
    }
    
    /** Set the detector size */
    void set_size(scitbx::vec2 <int> size) {
        size_ = size;
    }
    
    /* Set the distance from the crystal */
    void set_distance(double distance) {
        distance_ = distance;
    }
        
private:

    scitbx::vec3 <double> x_axis_;
    scitbx::vec3 <double> y_axis_;
    scitbx::vec3 <double> normal_;
    scitbx::vec2 <double> origin_;
    scitbx::vec2 <double> pixel_size_;
    scitbx::vec2 <int> size_;
    double distance_;
};

}} // namespace = dials::equipment

#endif // DIALS_EQUIPMENT_DETECTOR_H
