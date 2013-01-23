
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/format.hpp>
#include <string>
#include "../goniometer.h"

using namespace boost::python;
using namespace dials::equipment;

namespace dials { namespace equipment { 

namespace boost_python {

std::string goniometer_to_string(const Goniometer &goniometer) {
    boost::format fmt(
        "Goniometer:\n"
        "    rotation axis:     (%1%, %2%, %3%)\n"
        "    starting angle:    %4%\n"
        "    oscillation range: %5%\n"
        "    starting frame:    %6%\n"
        "    num frames:        %7%");
        
    fmt % goniometer.get_rotation_axis()[0];
    fmt % goniometer.get_rotation_axis()[1];
    fmt % goniometer.get_rotation_axis()[2];
    fmt % goniometer.get_starting_angle();
    fmt % goniometer.get_oscillation_range();
    fmt % goniometer.get_starting_frame();
    fmt % goniometer.get_num_frames();
    return fmt.str();
}

void export_goniometer() 
{
    class_ <Goniometer> ("Goniometer")
        .def(init <scitbx::vec3 <double>, double, double, int, int> ((
                arg("rotation_axis"), 
                arg("starting_angle"), 
                arg("oscillation_range"),
                arg("starting_frame"),
                arg("num_frames") = 0)))
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
        .add_property("num_frames",
            &Goniometer::get_num_frames,
            &Goniometer::set_num_frames)
        .def("get_angle_from_frame", 
            &Goniometer::get_angle_from_frame, (
                arg("frame")))
        .def("get_frame_from_angle",
            &Goniometer::get_frame_from_angle, (
                arg("angle")))
        .def("__str__", &goniometer_to_string);                
}

} // namespace = boost_python

}} // namespace = dials::equipment
