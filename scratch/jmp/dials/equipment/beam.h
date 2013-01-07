
#ifndef DIALS_EQUIPMENT_BEAM_H
#define DIALS_EQUIPMENT_BEAM_H

#include <exception>
#include <scitbx/vec3.h>

namespace dials {
namespace equipment {

class beam {

public:

    beam() : _wavelength(0.0) {}

    beam(scitbx::vec3 <double> direction) 
        : _wavelength(direction.length()),
          _direction(direction) {}
    
    beam(scitbx::vec3 <double> direction, double wavelength)
        : _wavelength(wavelength),
          _direction(direction) {  
        _direction.normalize();
        _direction *= _wavelength;
    }
    
    scitbx::vec3 <double> get_direction() {
        return _direction;
    }
    
    double get_wavelength() {
        return _wavelength;
    }
      
private:

    double _wavelength;
    scitbx::vec3 <double> _direction;
};

}
}

#endif // DIALS_EQUIPMENT_BEAM_H
