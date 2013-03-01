
#ifndef DIALS_ALGORITHMS_INTEGRATION_FROM_BEAM_VECTOR_TO_XDS_H
#define DIALS_ALGORITHMS_INTEGRATION_FROM_BEAM_VECTOR_TO_XDS_H

#include <scitbx/constants.h>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include "xds_coordinate_system.h"

namespace dials { namespace algorithms {

  using scitbx::constants::pi_180;
  using scitbx::vec3;

  /** Constant for scaling values */
//  static const double r2d = 1.0 / pi_180;
  static const double r2d = 1.0;

  /**
   * Class to represent a geometry transform from beam vector to XDS coordinates
   */
  class FromBeamVectorToXds {

  public:

    /**
     * Initialise the transform using the XDS coordinate system.
     * @param xcs The XDS coordinate system
     * @param s1 The diffracted beam vector
     * @param phi The rotation angle
     */
    FromBeamVectorToXds(XdsCoordinateSystem xcs,
                        vec3 <double> s1,
                        double phi)
      : scaled_e1_(xcs.get_e1_axis() * r2d / s1.length()),
        scaled_e2_(xcs.get_e2_axis() * r2d / s1.length()),
        s1_(s1),
        phi_(phi),
        zeta_(xcs.get_zeta()) {}

  public:

    /**
     * Apply the transform to a beam vector
     * @param s_dash The diffracted beam vector to transform
     * @param phi_dash The rotation angle for the beam vector.
     * @returns The point in XDS coordinates
     */
    vec3 <double> apply(vec3 <double> s_dash, double phi_dash) const {
      return vec3 <double> (
          scaled_e1_ * (s_dash - s1_),
          scaled_e2_ * (s_dash - s1_),
          zeta_ * (phi_dash - phi_));
    }

  private:

      vec3 <double> scaled_e1_;
      vec3 <double> scaled_e2_;
      vec3 <double> s1_;
      double phi_;
      double zeta_;
  };

}} // namespace = dials::algorithms

#endif // DIALS_ALGORITHMS_INTEGRATION_FROM_BEAM_VECTOR_TO_XDS_H
