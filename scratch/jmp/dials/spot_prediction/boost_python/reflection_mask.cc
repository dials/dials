
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../reflection_mask.h"

using namespace boost::python;

namespace dials { namespace spot_prediction { 

namespace boost_python {

void export_reflection_mask()
{
    class_ <ReflectionMask> ("ReflectionMask")
        .def(init <scitbx::vec3 <int>,
                   scitbx::vec3 <int> > ((
            arg("size"), 
            arg("roi_size"))))
        .def("create",
            &ReflectionMask::create, (
                arg("reflection_xyz")))            
        .add_property("mask", 
            &ReflectionMask::get_mask);
}

} // namespace = boost_python

}} // namespace = dials::spot_prediction
