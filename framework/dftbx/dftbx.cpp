#include <dftbx.h>

BOOST_PYTHON_MODULE(dftbx_ext)
{
  dftbx::init_module();
  boost::python::class_<dftbx::int_ol>("int_ol")
    .def("size", & dftbx::int_ol::size)
    .def("push_back", & dftbx::int_ol::push_back);
}
