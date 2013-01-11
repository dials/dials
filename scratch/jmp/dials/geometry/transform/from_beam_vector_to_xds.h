
#ifndef DIALS_GEOMETRY_TRANSFORM_FROM_BEAM_VECTOR_TO_XDS_H
#define DIALS_GEOMETRY_TRANSFORM_FROM_BEAM_VECTOR_TO_XDS_H

#include <scitbx/constants.h>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include "../xds_coordinate_system.h"

namespace dials { namespace geometry { namespace transform {

/** Constant for scaling values */
static const double r2d = 1.0 / scitbx::constants::pi_180;

/** 
 * Class to represent a geometry transform from beam vector to XDS coordinates
 */
class FromBeamVectorToXds {

public:

    /** Default constructor */
    FromBeamVectorToXds()
        : _phi(0.0),
          _zeta(0.0) {}
    
    /**
     * Initialise the transform using the XDS coordinate system.
     * @param xcs The XDS coordinate system
     * @param s1 The diffracted beam vector
     * @param phi The rotation angle
     */
    FromBeamVectorToXds(XdsCoordinateSystem xcs,
                        scitbx::vec3 <double> s1,
                        double phi)
        : _scaled_e1(xcs.get_e1_axis() * r2d / s1.length()),
          _scaled_e2(xcs.get_e2_axis() * r2d / s1.length()),
          _s1(s1),
          _phi(phi),
          _zeta(xcs.get_zeta()) {}

public:

    /**
     * Apply the transform to a beam vector
     * @param s_dash The diffracted beam vector to transform
     * @param phi_dash The rotation angle for the beam vector.
     * @returns The point in XDS coordinates
     */
    scitbx::vec3 <double> apply(scitbx::vec3 <double> s_dash, double phi_dash) {
        return scitbx::vec3 <double> (
            _scaled_e1 * (s_dash - _s1), 
            _scaled_e2 * (s_dash - _s1), 
            _zeta * (phi_dash - _phi));
    }

private:

    scitbx::vec3 <double> _scaled_e1;
    scitbx::vec3 <double> _scaled_e2;
    scitbx::vec3 <double> _s1;
    double _phi;
    double _zeta;
};

}}} // namespace = dials::geometry::transform

#endif // DIALS_GEOMETRY_TRANSFORM_FROM_BEAM_VECTOR_TO_XDS_H
