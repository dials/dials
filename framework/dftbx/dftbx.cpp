#include <dftbx.h>

BOOST_PYTHON_MODULE(dftbx_ext)
{
  dftbx::init_module();
  boost::python::class_<dftbx::int_ol>("int_ol")
    .def("size", & dftbx::int_ol::size)
    .def("push_back", & dftbx::int_ol::push_back)
    .def("duplicate", & dftbx::int_ol::duplicate)
    .def("merge", & dftbx::int_ol::merge);
  boost::python::class_<dftbx::int_rl>("int_rl")
    .def("merge", & dftbx::int_rl::merge);
}

