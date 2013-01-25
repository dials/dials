#include <x2tbx.h>

BOOST_PYTHON_MODULE(x2tbx_ext)
{
  x2tbx::init_module();
  boost::python::class_<x2tbx::ReflectionList>("ReflectionList")
    .def("setup", & x2tbx::ReflectionList::setup)
    .def("set_unit_cell", & x2tbx::ReflectionList::set_unit_cell)
    .def("get_indices", & x2tbx::ReflectionList::get_indices)
    .def("setup_resolution_shells",
         & x2tbx::ReflectionList::setup_resolution_shells)
    .def("get_shell", & x2tbx::ReflectionList::get_shell)
    .def("shell_high_limits", & x2tbx::ReflectionList::shell_high_limits)
    .def("shell_low_limits", & x2tbx::ReflectionList::shell_low_limits)
    .def("merge", & x2tbx::ReflectionList::merge)
    .def("i_sigma", & x2tbx::ReflectionList::i_sigma)
    .def("total_i_sigma", & x2tbx::ReflectionList::total_i_sigma)
    .def("rmerge", & x2tbx::ReflectionList::rmerge)
    .def("i_sigma_shells", & x2tbx::ReflectionList::i_sigma_shells)
    .def("total_i_sigma_shells", & x2tbx::ReflectionList::total_i_sigma_shells)
    .def("rmerge_shells", & x2tbx::ReflectionList::rmerge_shells);
}
