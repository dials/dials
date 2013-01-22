
#include <boost/python.hpp>
#include <boost/python/def.hpp>

using namespace boost::python;

namespace dials { namespace array_family { 

namespace boost_python {

void export_sorting();

BOOST_PYTHON_MODULE(dials_array_family_ext)
{
    export_sorting();
}

} // namespace = boost_python

}} // namespace = dials::array_family
