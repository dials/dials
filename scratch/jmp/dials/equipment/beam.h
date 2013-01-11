
#ifndef DIALS_EQUIPMENT_BEAM_H
#define DIALS_EQUIPMENT_BEAM_H

#include <exception>
#include <scitbx/vec3.h>

namespace dials { namespace equipment {

/** The beam */
class Beam {

public:

    /** Default constructor */
    Beam() : _wavelength(0.0) {}
    
    /** 
     * Initialise beam
     * @param direction The beam direction vector
     * @param wavelength The wavelength of the beam
     */    
    Beam(scitbx::vec3 <double> direction, double wavelength)
        : _wavelength(wavelength),
          _direction(direction) {}
    
    /** Get the beam direction */
    scitbx::vec3 <double> get_direction() {
        return _direction;
    }
    
    /** Get the beam wavelength */
    double get_wavelength() {
        return _wavelength;
    }
    
    /** Set the beam direction */
    void set_direction(scitbx::vec3 <double> direction) {
        _direction = direction;
    }
    
    /** Set the beam wavelength */
    void set_wavelength(double wavelength) {
        _wavelength = wavelength;
    }
      
private:

    double _wavelength;
    scitbx::vec3 <double> _direction;
};

}} // namespace = dials::equipment

#endif // DIALS_EQUIPMENT_BEAM_H
