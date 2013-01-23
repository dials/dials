
#include <boost/python.hpp>
#include <boost/python/def.hpp>

using namespace boost::python;

namespace dials { namespace af { 

namespace boost_python {

void export_sorting();
void export_remove();

BOOST_PYTHON_MODULE(dials_array_family_ext)
{
    export_sorting();
    export_remove();
}

} // namespace = boost_python

}} // namespace = dials::array_family
