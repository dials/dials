
#include <boost/python.hpp>
#include <boost/python/def.hpp>

using namespace boost::python;

namespace dials { namespace integration { 

namespace boost_python {

void export_reflection_mask();

BOOST_PYTHON_MODULE(dials_integration_ext)
{
    export_reflection_mask();
}

} // namespace = boost_python

}} // namespace = dials::integration
