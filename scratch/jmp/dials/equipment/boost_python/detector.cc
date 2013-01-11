
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../detector.h"

using namespace boost::python;
using namespace dials::equipment;

namespace dials { namespace equipment { 

namespace boost_python {

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
            &Detector::set_distance);
}

} // namespace = boost_python

}} // namespace = dials::equipment
