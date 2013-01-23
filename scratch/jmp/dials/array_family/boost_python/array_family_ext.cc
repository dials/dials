
#include <boost/python.hpp>
#include <boost/python/def.hpp>

using namespace boost::python;

namespace dials { namespace af { 

namespace boost_python {

void export_sorting();
void export_remove();
void export_flex_tiny();
void export_test_tiny_int();

BOOST_PYTHON_MODULE(dials_array_family_ext)
{
    export_flex_tiny();
    export_sorting();
    export_remove();
    export_test_tiny_int();
}

} // namespace = boost_python

}} // namespace = dials::array_family
