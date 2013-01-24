#include <x2tbx.h>

BOOST_PYTHON_MODULE(x2tbx_ext)
{
  x2tbx::init_module();
  boost::python::class_<x2tbx::ReflectionList>("ReflectionList")
    .def("setup", & x2tbx::ReflectionList::setup)
    .def("get_indices", & x2tbx::ReflectionList::get_indices)
    .def("setup_resolution_shells",
         & x2tbx::ReflectionList::setup_resolution_shells)
    .def("merge", & x2tbx::ReflectionList::merge)
    .def("i_sigma", & x2tbx::ReflectionList::i_sigma)
    .def("rmerge", & x2tbx::ReflectionList::rmerge);
}
