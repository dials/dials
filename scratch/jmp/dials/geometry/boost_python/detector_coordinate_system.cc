
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../detector_coordinate_system.h"

using namespace boost::python;

namespace dials { namespace geometry { 

namespace boost_python {

void export_detector_coordinate_system()
{
    class_ <detector_coordinate_system> ("detector_coordinate_system")
        .def(init <scitbx::vec3 <double>,
                   scitbx::vec3 <double> > ((
                arg("x_axis"), 
                arg("y_axis"))))
        .def(init <scitbx::vec3 <double>,
                   scitbx::vec3 <double>,
                   scitbx::vec3 <double> > ((
                arg("x_axis"), 
                arg("y_axis"), 
                arg("normal"))))
        .add_property("x_axis",  
            &detector_coordinate_system::get_x_axis,
            &detector_coordinate_system::set_x_axis)
        .add_property("y_axis", 
            &detector_coordinate_system::get_y_axis,
            &detector_coordinate_system::set_y_axis)
        .add_property("normal", 
            &detector_coordinate_system::get_normal,
            &detector_coordinate_system::set_normal)
        .def("in_pixel_units",
            &detector_coordinate_system::in_pixel_units, (
                arg("pixel_size")))
        .def("in_si_units",
            &detector_coordinate_system::in_si_units, (
                arg("pixel_size")));     
}

}

}}
