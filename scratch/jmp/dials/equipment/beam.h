
#ifndef DIALS_EQUIPMENT_BEAM_H
#define DIALS_EQUIPMENT_BEAM_H

#include <exception>
#include <scitbx/vec3.h>

namespace dials { namespace equipment {

/** The beam */
class Beam {

public:

    /** Default constructor */
    Beam() : wavelength_(0.0) {}
    
    /** 
     * Initialise beam
     * @param direction The beam direction vector
     * @param wavelength The wavelength of the beam
     */    
    Beam(scitbx::vec3 <double> direction, double wavelength)
        : wavelength_(wavelength),
          direction_(direction) {}
    
    /** Get the beam direction */
    scitbx::vec3 <double> get_direction() {
        return direction_;
    }
    
    /** Get the beam wavelength */
    double get_wavelength() {
        return wavelength_;
    }
    
    /** Set the beam direction */
    void set_direction(scitbx::vec3 <double> direction) {
        direction_ = direction;
    }
    
    /** Set the beam wavelength */
    void set_wavelength(double wavelength) {
        wavelength_ = wavelength;
    }
      
private:

    double wavelength_;
    scitbx::vec3 <double> direction_;
};

}} // namespace = dials::equipment

#endif // DIALS_EQUIPMENT_BEAM_H
