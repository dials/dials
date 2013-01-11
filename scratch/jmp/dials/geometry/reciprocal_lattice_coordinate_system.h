
#ifndef DIALS_GEOMETRY_RECIPROCAL_LATTICE_COORDINATE_SYSTEM_H
#define DIALS_GEOMETRY_RECIPROCAL_LATTICE_COORDINATE_SYSTEM_H

#include <scitbx/vec3.h>
#include <scitbx/mat3.h>

namespace dials { namespace geometry {

/** Class representing reciprocal lattice coordinate system */
class ReciprocalLatticeCoordinateSystem {

public:

    /** Default constructor */
    ReciprocalLatticeCoordinateSystem() {}

    /** 
     * Initialise coordinate system using axis vectors
     * @param b1_star The h component vector
     * @param b2_star The k component vector
     * @param b3_star The l component vector
     */
    ReciprocalLatticeCoordinateSystem(scitbx::vec3 <double> b1_star,
                                      scitbx::vec3 <double> b2_star,
                                      scitbx::vec3 <double> b3_star)
        : b1_star_(b1_star),
          b2_star_(b2_star),
          b3_star_(b3_star) {}
    
    /** 
     * Initialise the coordinate system from a UB matrix
     * @param ub The UB matrix
     */
    ReciprocalLatticeCoordinateSystem(scitbx::mat3 <double> ub)
        : b1_star_(ub[0], ub[3], ub[6]),
          b2_star_(ub[1], ub[4], ub[7]),
          b3_star_(ub[2], ub[5], ub[8]) {}

public:

    /** Get the b1* axis vector */
    scitbx::vec3 <double> get_b1_star_axis() {
        return b1_star_;
    }
    
    /** Get the b2* axis vector */
    scitbx::vec3 <double> get_b2_star_axis() {
        return b2_star_;
    }
    
    /** Get the b3* axis vector */
    scitbx::vec3 <double> get_b3_star_axis() {
        return b3_star_;
    }
    
    /** Set the b1* axis vector */
    void set_b1_star_axis(scitbx::vec3 <double> b1_star) {
        b1_star_ = b1_star;
    }
    
    /** Set the b2* axis vector */
    void set_b2_star_axis(scitbx::vec3 <double> b2_star) {
        b2_star_ = b2_star;
    }
    
    /** Set the b3* axis vector */
    void set_b3_star_axis(scitbx::vec3 <double> b3_star) {
        b3_star_ = b3_star;
    }
    
    /**
     * Read the coordinate system from the UB matrix
     * @param ub The ub matrix
     */
    void from_ub_matrix(scitbx::mat3 <double> ub) {
        b1_star_ = scitbx::vec3 <double> (ub[0], ub[3], ub[6]);
        b2_star_ = scitbx::vec3 <double> (ub[1], ub[4], ub[7]);
        b3_star_ = scitbx::vec3 <double> (ub[2], ub[5], ub[8]);
    }
    
    /**
     * Convert the coordinate system to a UB matrix
     * @returns The UB matrix
     */
    scitbx::mat3 <double> to_ub_matrix() {
        return scitbx::mat3 <double> (
            b1_star_[0], b2_star_[0], b3_star_[0],
            b1_star_[1], b2_star_[1], b3_star_[1],
            b1_star_[2], b2_star_[2], b3_star_[2]
        );
    }

private:

    scitbx::vec3 <double> b1_star_;
    scitbx::vec3 <double> b2_star_;
    scitbx::vec3 <double> b3_star_;
};

}} // namespace = dials::geometry

#endif // DIALS_GEOMETRY_RECIPROCAL_LATTICE_COORDINATE_SYSTEM_H
