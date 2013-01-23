
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/format.hpp>
#include <string>
#include "../beam.h"

using namespace boost::python;

namespace dials { namespace equipment { 

namespace boost_python {

std::string beam_to_string(const Beam &beam) {
    boost::format fmt(
        "Beam:\n"
        "    wavelength: %1%\n"
        "    direction : (%2%, %3%, %4%)");
        
    fmt % beam.get_wavelength();
    fmt % beam.get_direction()[0];
    fmt % beam.get_direction()[1];
    fmt % beam.get_direction()[2];
    return fmt.str();
}

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
            &Beam::set_wavelength)
        .def("__str__", &beam_to_string);
}

} // namespace = boost_python

}} // namespace = dials::equipment
