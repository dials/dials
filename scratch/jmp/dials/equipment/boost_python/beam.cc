
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../beam.h"

using namespace boost::python;

namespace dials { namespace equipment { 

namespace boost_python {

void export_beam()
{
    class_ <beam> ("beam")
        .def(init <scitbx::vec3 <double>, double> ((
            arg("direction"), 
            arg("wavelength"))))
        .add_property("direction",  
            &beam::get_direction,
            &beam::set_direction)
        .add_property("wavelength", 
            &beam::get_wavelength,
            &beam::set_wavelength);
}

}

}}
