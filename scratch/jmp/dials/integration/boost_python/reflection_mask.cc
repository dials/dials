
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../reflection_mask.h"

using namespace boost::python;

namespace dials { namespace integration { 

namespace boost_python {

void export_reflection_mask()
{
    scitbx::af::flex_bool (ReflectionMask::*create_arrays)(
                const dials::af::flex_vec3_double &,
                const dials::af::flex_tiny6_int &) = &ReflectionMask::create;

    scitbx::af::flex_bool (ReflectionMask::*create_list)(
                ReflectionList &) = &ReflectionMask::create;

    class_ <ReflectionMask> ("ReflectionMask")
        .def(init <scitbx::vec3 <int> > ((
            arg("mask_size"))))
        .def("create",
            create_arrays, (
                arg("xyz"),
                arg("roi")))
        .def("create",                
            create_list, (
                arg("reflections")))
        .add_property("mask", 
            &ReflectionMask::get_mask);
}

} // namespace = boost_python

}} // namespace = dials::integration
