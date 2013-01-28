
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/format.hpp>
#include <string>
#include "../detector.h"

using namespace boost::python;
using namespace dials::equipment;

namespace dials { namespace equipment { 

namespace boost_python {

std::string detector_to_string(const Detector &detector) {
    boost::format fmt(
        "Detector:\n"
        "    x axis:     (%1%, %2%, %3%)\n"
        "    y axis:     (%4%, %5%, %6%)\n"
        "    normal:     (%7%, %8%, %9%)\n"
        "    origin:     (%10%, %11%)\n"
        "    size:       (%12%, %13%)\n"
        "    pixel size: (%14%, %15%)\n"
        "    distance:   (%16%)");
        
    fmt % detector.get_x_axis()[0];
    fmt % detector.get_x_axis()[1];
    fmt % detector.get_x_axis()[2];
    fmt % detector.get_y_axis()[0];
    fmt % detector.get_y_axis()[1];
    fmt % detector.get_y_axis()[2];
    fmt % detector.get_normal()[0];
    fmt % detector.get_normal()[1];
    fmt % detector.get_normal()[2];
    fmt % detector.get_origin()[0];
    fmt % detector.get_origin()[1];
    fmt % detector.get_size()[0];
    fmt % detector.get_size()[1];
    fmt % detector.get_pixel_size()[0];
    fmt % detector.get_pixel_size()[1];  
    fmt % detector.get_distance();  
    return fmt.str();
}

void export_detector() 
{
    class_ <Detector> ("Detector")
        .def(init <scitbx::vec3 <double>, 
                   scitbx::vec3 <double>,
                   scitbx::vec3 <double>, 
                   scitbx::vec2 <double>,
                   scitbx::vec2 <double>,
                   scitbx::vec2 <int>, 
                   double> ((
                arg("x_axis"), 
                arg("y_axis"), 
                arg("normal"), 
                arg("origin"), 
                arg("pixel_size"),
                arg("size"),
                arg("distance"))))
        .add_property("x_axis",  
            &Detector::get_x_axis,
            &Detector::set_x_axis)
        .add_property("y_axis",  
            &Detector::get_y_axis,
            &Detector::set_y_axis)
        .add_property("normal",  
            &Detector::get_normal,
            &Detector::set_normal)
        .add_property("origin",  
            &Detector::get_origin,
            &Detector::set_origin)
        .add_property("pixel_size",  
            &Detector::get_pixel_size,
            &Detector::set_pixel_size)
        .add_property("size",
            &Detector::get_size,
            &Detector::set_size)
        .add_property("distance",  
            &Detector::get_distance,
            &Detector::set_distance)
        .def("is_coordinate_valid",
            &Detector::is_coordinate_valid, (
                arg("xy")))
        .def("__str__", &detector_to_string);            
}

} // namespace = boost_python

}} // namespace = dials::equipment
