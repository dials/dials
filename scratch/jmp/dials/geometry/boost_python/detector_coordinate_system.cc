
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../detector_coordinate_system.h"

using namespace boost::python;
using namespace dials::geometry;

namespace dials { namespace geometry { namespace boost_python {

void export_detector_coordinate_system()
{
    class_<detector_coordinate_system>("detector_coordinate_system")
        .def(init <scitbx::vec3 <double>,
                   scitbx::vec3 <double>,
                   scitbx::vec3 <double> > ((
            arg("axis_x"), arg("axis_y"), arg("normal"))))
        .def(init <scitbx::vec3 <double>,
                   scitbx::vec3 <double> > ((
            arg("axis_x"), arg("axis_y"))))
        .add_property("axis_x",  
            &detector_coordinate_system::get_axis_x,
            &detector_coordinate_system::set_axis_x)
        .add_property("axis_y", 
            &detector_coordinate_system::get_axis_y,
            &detector_coordinate_system::set_axis_y)
        .add_property("normal", 
            &detector_coordinate_system::get_normal,
            &detector_coordinate_system::set_normal)
        .def("in_pixel_units",
            &detector_coordinate_system::in_pixel_units,
            (arg("pixel_size")))
        .def("in_si_units",
            &detector_coordinate_system::in_si_units,
            (arg("pixel_size")))
        .staticmethod("in_pixel_units")
        .staticmethod("in_si_units");
}

}}}
