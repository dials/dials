
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../reflection_mask.h"

using namespace boost::python;

namespace dials { namespace integration { 

namespace boost_python {

void export_reflection_mask()
{
    class_ <ReflectionMask> ("ReflectionMask")
        .def(init <scitbx::vec3 <int> > ((
            arg("mask_size"))))
        .def("create",
            &ReflectionMask::create, (
                arg("xyz"),
                arg("roi")))
        .add_property("mask", 
            &ReflectionMask::get_mask);
}

} // namespace = boost_python

}} // namespace = dials::integration
