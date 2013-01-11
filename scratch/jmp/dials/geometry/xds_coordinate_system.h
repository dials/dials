
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
        : e1_(s1.cross(s0).normalize()),
          e2_(s1.cross(e1_).normalize()),
          e3_((s1 + s0).normalize()),
          zeta_(m2 * e1_) {}

public:

    /** Get the e1 axis vector */
    scitbx::vec3 <double> get_e1_axis() {
        return e1_;
    }
    
    /** Get the e2 axis vector */
    scitbx::vec3 <double> get_e2_axis() {
        return e2_;
    }
    
    /** Get the e3 axis vector */
    scitbx::vec3 <double> get_e3_axis() {
        return e3_;
    }
    
    /** Get the lorentz correction factor (zeta) */
    double get_zeta() {
        return zeta_;
    }

private:

    scitbx::vec3 <double> e1_;
    scitbx::vec3 <double> e2_;
    scitbx::vec3 <double> e3_;
    double zeta_;
};

}} // namespace = dials::geometry

#endif // DIALS_GEOMETRY_XDS_COORDINATE_SYSTEM_H
