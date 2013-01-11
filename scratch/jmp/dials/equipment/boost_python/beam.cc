
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../beam.h"

using namespace boost::python;

namespace dials { namespace equipment { 

namespace boost_python {

void export_beam()
{
    class_ <Beam> ("Beam")
        .def(init <scitbx::vec3 <double>, double> ((
            arg("direction"), 
            arg("wavelength"))))
        .add_property("direction",  
            &Beam::get_direction,
            &Beam::set_direction)
        .add_property("wavelength", 
            &Beam::get_wavelength,
            &Beam::set_wavelength);
}

} // namespace = boost_python

}} // namespace = dials::equipment
