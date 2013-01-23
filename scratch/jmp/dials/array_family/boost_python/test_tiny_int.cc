
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <scitbx/array_family/tiny.h>

using namespace boost::python;

namespace dials { namespace af {

namespace boost_python {

scitbx::af::tiny <int, 3> test_tiny_int3() {
    return scitbx::af::tiny <int, 3> (
        1, 2, 3);
}

scitbx::af::tiny <int, 4> test_tiny_int4() {
    return scitbx::af::tiny <int, 4> (
        1, 2, 3, 4);
}

scitbx::af::tiny <int, 6> test_tiny_int6() {
    return scitbx::af::tiny <int, 6> (
        1, 2, 3, 4, 5, 6);
}

scitbx::af::tiny <double, 6> test_tiny_double6() {
    return scitbx::af::tiny <double, 6> (
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
}

void export_test_tiny_int() 
{
    def("test_tiny_int3", &test_tiny_int3);
    def("test_tiny_int4", &test_tiny_int4);
    def("test_tiny_int6", &test_tiny_int6);
    def("test_tiny_double6", &test_tiny_double6); 
}

}

}}
