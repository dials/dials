
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../detector.h"

using namespace boost::python;
using namespace dials::equipment;

namespace dials { namespace equipment { namespace boost_python {

void export_detector() 
{
    class_<detector>("detector")
        .def(init <scitbx::vec3 <double>, scitbx::vec3 <double>,
                   scitbx::vec3 <double>, scitbx::vec2 <double>,
                   scitbx::vec2 <double>, double> ((
            arg("axis_x"), arg("axis_y"), arg("normal"), 
            arg("origin"), arg("pixel_size"), arg("distance"))))
        .add_property("axis_x",  
            &detector::get_axis_x,
            &detector::set_axis_x)
        .add_property("axis_y",  
            &detector::get_axis_y,
            &detector::set_axis_y)
        .add_property("normal",  
            &detector::get_normal,
            &detector::set_normal)
        .add_property("origin",  
            &detector::get_origin,
            &detector::set_origin)
        .add_property("pixel_size",  
            &detector::get_pixel_size,
            &detector::set_pixel_size)
        .add_property("distance",  
            &detector::get_distance,
            &detector::set_distance);
}

}}}
