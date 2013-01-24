#include <x2tbx.h>

BOOST_PYTHON_MODULE(x2tbx_ext)
{
  x2tbx::init_module();
  boost::python::class_<x2tbx::resolutionizer>("resolutionizer")
    .def("set_unit_cell", & x2tbx::resolutionizer::set_unit_cell)
    .def("compare_resolution", & x2tbx::resolutionizer::compare_resolution)
    .def("setup", & x2tbx::resolutionizer::setup)
    .def("sorted_indices", & x2tbx::resolutionizer::sorted_indices)
    .def("setup_shells", & x2tbx::resolutionizer::setup_shells)
    .def("isig_shells", & x2tbx::resolutionizer::isig_shells);
  boost::python::class_<x2tbx::ReflectionList>("ReflectionList")
    .def("setup", & x2tbx::ReflectionList::setup)
    .def("merge", & x2tbx::ReflectionList::merge)
    .def("i_sigma", & x2tbx::ReflectionList::i_sigma)
    .def("rmerge", & x2tbx::ReflectionList::rmerge);
}
