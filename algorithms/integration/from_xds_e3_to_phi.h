/*
 * from_xds_e3_to_phi.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_INTEGRATION_FROM_XDS_E3_TO_PHI_H
#define DIALS_ALGORITHMS_INTEGRATION_FROM_XDS_E3_TO_PHI_H

namespace dials { namespace algorithms {

  /** A class to transform from XDS e3 coord to the rotation angle, phi */
  class FromXdsE3ToPhi {

  public:

    /**
     * Initialise the class
     * @param zeta The xds zeta parameters
     * @param phi The rotation angle
     */
    FromXdsE3ToPhi(double zeta, double phi)
      : zeta_(zeta),
        phi_(phi) {}

    /**
     * Apply the transform
     * @param e3 The XDS e3 coordinate
     * @returns The rotation angle phi'
     */
    double operator()(double e3) const {
      return e3 / zeta_ + phi_;
    }

  private:

    double zeta_;
    double phi_;
  };

}} // namespace dials::algorithms

#endif // DIALS_ALGORITHMS_INTEGRATION_FROM_XDS_E3_TO_PHI_H
