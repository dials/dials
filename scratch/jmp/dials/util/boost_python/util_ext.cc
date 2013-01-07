
#include <exception>
#include <boost/python.hpp>

using namespace boost::python;

namespace dials { namespace util { namespace boost_python {

void std_exception_translator(std::exception const& x) {
    PyErr_SetString(PyExc_UserWarning, x.what());
}


BOOST_PYTHON_MODULE(util_ext)
{
    register_exception_translator<std::exception>(std_exception_translator);
}

}}}
