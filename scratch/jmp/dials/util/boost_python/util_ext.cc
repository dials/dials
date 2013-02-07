
#include <boost/python.hpp>
#include <exception>

using namespace boost::python;

namespace dials { namespace util { 

namespace boost_python {

void export_test_inheritance();

void std_exception_translator(std::exception const& x) {
    PyErr_SetString(PyExc_UserWarning, x.what());
}


BOOST_PYTHON_MODULE(dials_util_ext)
{
    register_exception_translator<std::exception>(std_exception_translator);
    export_test_inheritance();
}

} // namespace = boost_python

}} // namespace = dials::util
