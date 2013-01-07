
#ifndef DIALS_GEOMETRY_RECIPROCAL_LATTICE_COORDINATE_SYSTEM_H
#define DIALS_GEOMETRY_RECIPROCAL_LATTICE_COORDINATE_SYSTEM_H

#include <scitbx/vec3.h>
#include <scitbx/mat3.h>

namespace dials { namespace geometry {

/** Class representing reciprocal lattice coordinate system */
class reciprocal_lattice_coordinate_system {

public:

    /** Default constructor */
    reciprocal_lattice_coordinate_system() {}

    /** 
     * Initialise coordinate system using axis vectors
     * @param b1_star The h component vector
     * @param b2_star The k component vector
     * @param b3_star The l component vector
     */
    reciprocal_lattice_coordinate_system(scitbx::vec3 <double> b1_star,
                                         scitbx::vec3 <double> b2_star,
                                         scitbx::vec3 <double> b3_star)
        : _b1_star(b1_star),
          _b2_star(b2_star),
          _b3_star(b3_star) {}
    
    /** 
     * Initialise the coordinate system from a UB matrix
     * @param ub The UB matrix
     */
    reciprocal_lattice_coordinate_system(scitbx::mat3 <double> ub)
        : _b1_star(ub[0], ub[1], ub[2]),
          _b2_star(ub[3], ub[4], ub[5]),
          _b3_star(ub[6], ub[7], ub[8]) {}

public:

    /** Get the b1* axis vector */
    scitbx::vec3 <double> get_b1_star() {
        return _b1_star;
    }
    
    /** Get the b2* axis vector */
    scitbx::vec3 <double> get_b2_star() {
        return _b2_star;
    }
    
    /** Get the b3* axis vector */
    scitbx::vec3 <double> get_b3_star() {
        return _b3_star;
    }
    
    /** Set the b1* axis vector */
    void set_b1_star(scitbx::vec3 <double> b1_star) {
        _b1_star = b1_star;
    }
    
    /** Set the b2* axis vector */
    void set_b2_star(scitbx::vec3 <double> b2_star) {
        _b2_star = b2_star;
    }
    
    /** Set the b3* axis vector */
    void set_b3_star(scitbx::vec3 <double> b3_star) {
        _b3_star = b3_star;
    }

private:

    scitbx::vec3 <double> _b1_star;
    scitbx::vec3 <double> _b2_star;
    scitbx::vec3 <double> _b3_star;
};

}}

#endif // DIALS_GEOMETRY_RECIPROCAL_LATTICE_COORDINATE_SYSTEM_H
