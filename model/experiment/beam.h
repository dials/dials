/*
 * beam.h
 *
 *   Copyright (C) 2013 Diamond Light Source, James Parkhurst
 *
 *   This code is distributed under the BSD license, a copy of which is
 *   included in the root directory of this package.
 */
#ifndef DIALS_MODEL_EXPERIMENT_BEAM_H
#define DIALS_MODEL_EXPERIMENT_BEAM_H

#include <scitbx/vec3.h>

namespace dials { namespace model { namespace experiment {

  using scitbx::vec3;

  /**
   * A class to represent the X-ray primary beam for a standard rotation
   * geometry diffraction experiment. The following assumptions are made:
   * (i)   the beam is monochromatic
   * (ii)  the beam is reasonably parallel
   * (iii) the direction of the beam vector is towards the sample
   * (iv)  the length of the direction vector is 1.0 / wavelength.
   *
   * In the document detailing the conventions used:
   *    direction -> s0
   *    wavelength -> lambda
   *    polarization -> *unspecified*
   *    polarization_fraction -> *unspecified*
   */
  class Beam {
  public:
    /** Default constructor: initialise all to zero */
    Beam()
      : direction_(0.0, 0.0, 0.0),
        wavelength_(0.0),
        polarization_(0.0, 0.0, 0.0),
        polarization_fraction_(0.0) {}

    /**
     * Initialise the beam object on the principle that the beam is aligned
     * with the +z axis, as is quite normal. Also assume the beam has
     * polarization fraction 0.999 and is polarized in the x-z plane.
     * @param wavelength The wavelength of the beam
     */
    Beam(double wavelength)
      : direction_(vec3 <double> (0.0, 0.0, 1.0) / wavelength),
        wavelength_(wavelength),
        polarization_(0.0, 1.0, 0.0),
        polarization_fraction_(0.999) {}

    /**
     * Initialise a beam object with a direction. Normalize the direction
     * vector and give it the length of 1.0 / wavelength. Assume the beam has
     * polarization fraction 0.999 and is polarized in the x-z plane.
     * @param wavelength The wavelength of the beam
     * @param direction The beam direction vector.
     */
    Beam(double wavelength,
         vec3 <double> direction)
      : direction_(direction / wavelength),
        wavelength_(wavelength),
        polarization_(0.0, 1.0, 0.0),
        polarization_fraction_(0.999) {}

    /**
     * Initialise all the beam parameters. Normalize the direction vector
     * and give it the length of 1.0 / wavelength
     * @param wavelength The wavelength of the beam
     * @param direction The beam direction vector.
     * @param polarization The polarization plane of the beam
     * @param polarization_fraction The polarization fraction.
     */
    Beam(double wavelength,
         vec3 <double> direction,
         vec3 <double> polarization,
         double polarization_fraction)
      : direction_(direction),
        wavelength_(1.0 / direction.length()),
        polarization_(polarization),
        polarization_fraction_(polarization_fraction) {}

    /**
     * Initialise a beam object with a direction. Assign the wavelength of the
     * beam to be 1.0 / length of the direction vector. Assume the beam has
     * polarization fraction 0.999 and is polarized in the x-z plane.
     * @param wavelength The wavelength of the beam
     * @param direction The beam direction vector.
     */
    Beam(vec3 <double> direction)
      : direction_(direction),
        wavelength_(1.0 / direction.length()),
        polarization_(0.0, 1.0, 0.0),
        polarization_fraction_(0.999) {}

    /**
     * Initialise all the beam parameters.
     * @param direction The beam direction vector.
     * @param polarization The polarization plane of the beam
     * @param polarization_fraction The polarization fraction.
     */
    Beam(vec3 <double> direction,
         vec3 <double> polarization,
         double polarization_fraction)
      : direction_(direction),
        wavelength_(1.0 / direction.length()),
        polarization_(polarization),
        polarization_fraction_(polarization_fraction) {}

    /** Get the direction */
    vec3 <double> get_direction() const {
      return direction_;
    }

    /** Get the wavelength */
    double get_wavelength() const {
      return wavelength_;
    }

    /** Get the polarization */
    vec3 <double> get_polarization() const {
      return polarization_;
    }

    /** Get the polarization fraction */
    double get_polarization_fraction() const {
      return polarization_fraction_;
    }

    /** Set the direction. */
    void set_direction(vec3 <double> direction) {
      direction_ = direction;
      wavelength_ = 1.0 / direction.length();
    }

    /** Set the polarization plane */
    void set_polarization(vec3 <double> polarization) {
      polarization_ = polarization;
    }

    /** Set the polarization fraction */
    void set_polarization_fraction(double polarization_fraction) {
      polarization_fraction_ = polarization_fraction;
    }

  private:
    vec3 <double> direction_;
    double wavelength_;
    vec3 <double> polarization_;
    double polarization_fraction_;
};

}}} // namespace dials::model::experiment

#endif // DIALS_MODEL_EXPERIMENT_BEAM_H
