
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../goniometer.h"

using namespace boost::python;
using namespace dials::equipment;

namespace dials { namespace equipment { 

namespace boost_python {

void export_goniometer() 
{
    class_ <Goniometer> ("Goniometer")
        .def(init <scitbx::vec3 <double>, double, double, int> ((
                arg("rotation_axis"), 
                arg("starting_angle"), 
                arg("oscillation_range"),
                arg("starting_frame"))))
        .add_property("rotation_axis",  
            &Goniometer::get_rotation_axis,
            &Goniometer::set_rotation_axis)
        .add_property("starting_angle",  
            &Goniometer::get_starting_angle,
            &Goniometer::set_starting_angle)
        .add_property("oscillation_range",  
            &Goniometer::get_oscillation_range,
            &Goniometer::set_oscillation_range)
        .add_property("starting_frame",  
            &Goniometer::get_starting_frame,
            &Goniometer::set_starting_frame)
        .def("get_angle_from_frame", 
            &Goniometer::get_angle_from_frame, (
                arg("frame")))
        .def("get_frame_from_angle",
            &Goniometer::get_frame_from_angle, (
                arg("angle")));
}

} // namespace = boost_python

}} // namespace = dials::equipment
