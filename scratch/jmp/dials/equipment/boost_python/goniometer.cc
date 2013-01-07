
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../goniometer.h"

using namespace boost::python;
using namespace dials::equipment;

namespace dials { namespace equipment { namespace boost_python {

void export_goniometer() 
{
    class_<goniometer>("goniometer")
        .def(init <scitbx::vec3 <double>, int, double, double> ((
            arg("rotation_axis"), arg("starting_frame"), 
            arg("starting_angle"), arg("oscillation_range"))))
        .add_property("rotation_axis",  
            &goniometer::get_rotation_axis,
            &goniometer::set_rotation_axis)
        .add_property("starting_frame",  
            &goniometer::get_starting_frame,
            &goniometer::set_starting_frame)
        .add_property("starting_angle",  
            &goniometer::get_starting_angle,
            &goniometer::set_starting_angle)
        .add_property("oscillation_range",  
            &goniometer::get_oscillation_range,
            &goniometer::set_oscillation_range)
        .def("get_angle_from_frame", 
            &goniometer::get_angle_from_frame, (arg("frame")))
        .def("get_frame_from_angle",
            &goniometer::get_frame_from_angle, (arg("angle")));
}

}}}
