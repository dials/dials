
#include <boost/python.hpp>
#include <boost/python/def.hpp>

using namespace boost::python;

namespace dials {

namespace boost_python {

void export_reflection();

BOOST_PYTHON_MODULE(dials_reflection_ext)
{
    export_reflection();
}

} // namespace = boost_python

} // namespace = dials
