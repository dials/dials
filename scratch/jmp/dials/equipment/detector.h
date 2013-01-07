
#ifndef DIALS_EQUIPMENT_DETECTOR_H
#define DIALS_EQUIPMENT_DETECTOR_H

#include <scitbx/vec2.h>
#include <scitbx/vec3.h>

namespace dials { namespace equipment {

/** The detector */
class detector {

public:

    /** The default constructor */
    detector() : _distance(0.0) {}
    
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
    detector(scitbx::vec3 <double> x_axis,
             scitbx::vec3 <double> y_axis,
             scitbx::vec3 <double> normal,
             scitbx::vec2 <double> origin,
             scitbx::vec2 <double> pixel_size,
             scitbx::vec2 <int> size,
             double distance)
        : _x_axis(x_axis),
          _y_axis(y_axis),
          _normal(normal),
          _origin(origin),
          _pixel_size(pixel_size),
          _size(size),
          _distance(distance) {}

public:

    /** Get the x axis */
    scitbx::vec3 <double> get_x_axis() {
        return _x_axis;
    }
    
    /** Get the y axis */
    scitbx::vec3 <double> get_y_axis() {
        return _y_axis;
    }
    
    /** Get the normal */
    scitbx::vec3 <double> get_normal() {
        return _normal;
    }

    /** Get the detector origin */
    scitbx::vec2 <double> get_origin() {
        return _origin;
    }
    
    /** Get the pixel size */
    scitbx::vec2 <double> get_pixel_size() {
        return _pixel_size;
    }

    /** Get the detector size */
    scitbx::vec2 <int> get_size() {
        return _size;
    }
    
    /** Get the distance from the crystal */
    double get_distance() {
        return _distance;
    }
    
    /** Set the x axis */
    void set_x_axis(scitbx::vec3 <double> x_axis) {
        _x_axis = x_axis;
    }
    
    /** Set the y axis */
    void set_y_axis(scitbx::vec3 <double> y_axis) {
        _y_axis = y_axis;
    }
    
    /** Set the normal */
    void set_normal(scitbx::vec3 <double> normal) {
        _normal = normal;
    }
        
    /** Set the origin */
    void set_origin(scitbx::vec2 <double> origin) {
        _origin = origin;
    }
    
    /** Set the pixel size */
    void set_pixel_size(scitbx::vec2 <double> pixel_size) {
        _pixel_size = pixel_size;
    }
    
    /** Set the detector size */
    void set_size(scitbx::vec2 <int> size) {
        _size = size;
    }
    
    /* Set the distance from the crystal */
    void set_distance(double distance) {
        _distance = distance;
    }
        
private:

    scitbx::vec3 <double> _x_axis;
    scitbx::vec3 <double> _y_axis;
    scitbx::vec3 <double> _normal;
    scitbx::vec2 <double> _origin;
    scitbx::vec2 <double> _pixel_size;
    scitbx::vec2 <int> _size;
    double _distance;
};

}}

#endif // DIALS_EQUIPMENT_DETECTOR_H
