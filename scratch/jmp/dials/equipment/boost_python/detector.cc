
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../detector.h"

using namespace boost::python;
using namespace dials::equipment;

namespace dials { namespace equipment { 

namespace boost_python {

void export_detector() 
{
    class_ <detector> ("detector")
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
            &detector::get_x_axis,
            &detector::set_x_axis)
        .add_property("y_axis",  
            &detector::get_y_axis,
            &detector::set_y_axis)
        .add_property("normal",  
            &detector::get_normal,
            &detector::set_normal)
        .add_property("origin",  
            &detector::get_origin,
            &detector::set_origin)
        .add_property("pixel_size",  
            &detector::get_pixel_size,
            &detector::set_pixel_size)
        .add_property("size",
            &detector::get_size,
            &detector::set_size)
        .add_property("distance",  
            &detector::get_distance,
            &detector::set_distance);
}

}

}}
