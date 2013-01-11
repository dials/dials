
#ifndef DIALS_GEOMETRY_XDS_COORDINATE_SYSTEM_H
#define DIALS_GEOMETRY_XDS_COORDINATE_SYSTEM_H

#include <scitbx/vec3.h>
#include <scitbx/mat3.h>

namespace dials { namespace geometry {

/** Class representing XDS coordinate system */
class XdsCoordinateSystem {

public:

    /** Default constructor */
    XdsCoordinateSystem() {}

    /** 
     * Initialise coordinate system
     * @param s0 The incident beam vector
     * @param s1 The diffracted beam vector
     * @param m2 The rotation axis
     * @param phi The rotation angle
     */
    XdsCoordinateSystem(scitbx::vec3 <double> s0,
                        scitbx::vec3 <double> s1,
                        scitbx::vec3 <double> m2,
                        double phi)
        : _e1(s1.cross(s0).normalize()),
          _e2(s1.cross(_e1).normalize()),
          _e3((s1 + s0).normalize()),
          _zeta(m2 * _e1) {}

public:

    /** Get the e1 axis vector */
    scitbx::vec3 <double> get_e1_axis() {
        return _e1;
    }
    
    /** Get the e2 axis vector */
    scitbx::vec3 <double> get_e2_axis() {
        return _e2;
    }
    
    /** Get the e3 axis vector */
    scitbx::vec3 <double> get_e3_axis() {
        return _e3;
    }
    
    /** Get the lorentz correction factor (zeta) */
    double get_zeta() {
        return _zeta;
    }

private:

    scitbx::vec3 <double> _e1;
    scitbx::vec3 <double> _e2;
    scitbx::vec3 <double> _e3;
    double _zeta;
};

}} // namespace = dials::geometry

#endif // DIALS_GEOMETRY_XDS_COORDINATE_SYSTEM_H
