
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <cctbx/miller.h>
#include "../remove.h"
#include "../array_types.h"

using namespace boost::python;

namespace dials { namespace af { 

namespace boost_python {

void export_remove()
{
    def("remove_if", &remove_if <scitbx::af::flex_bool>);
    def("remove_if", &remove_if <scitbx::af::flex_int>);
    def("remove_if", &remove_if <scitbx::af::flex_long>);
    def("remove_if", &remove_if <scitbx::af::flex_size_t>);
    def("remove_if", &remove_if <scitbx::af::flex_float>);
    def("remove_if", &remove_if <scitbx::af::flex_double>);
    def("remove_if", &remove_if <scitbx::af::flex_complex_double>);
    def("remove_if", &remove_if <dials::af::flex_vec2_int>);
    def("remove_if", &remove_if <dials::af::flex_vec2_double>);
    def("remove_if", &remove_if <dials::af::flex_vec3_int>);
    def("remove_if", &remove_if <dials::af::flex_vec3_double>);
    def("remove_if", &remove_if <dials::af::flex_tiny6_int>);
    def("remove_if", &remove_if <dials::af::flex_tiny6_double>);
    def("remove_if", &remove_if <dials::af::flex_miller_index>);

    def("remove_if_not", &remove_if_not <scitbx::af::flex_bool>);
    def("remove_if_not", &remove_if_not <scitbx::af::flex_int>);
    def("remove_if_not", &remove_if_not <scitbx::af::flex_long>);
    def("remove_if_not", &remove_if_not <scitbx::af::flex_size_t>);
    def("remove_if_not", &remove_if_not <scitbx::af::flex_float>);
    def("remove_if_not", &remove_if_not <scitbx::af::flex_double>);
    def("remove_if_not", &remove_if_not <scitbx::af::flex_complex_double>);
    def("remove_if_not", &remove_if_not <dials::af::flex_vec2_int>);
    def("remove_if_not", &remove_if_not <dials::af::flex_vec2_double>);
    def("remove_if_not", &remove_if_not <dials::af::flex_vec3_int>);
    def("remove_if_not", &remove_if_not <dials::af::flex_vec3_double>);
    def("remove_if_not", &remove_if_not <dials::af::flex_tiny6_int>);
    def("remove_if_not", &remove_if_not <dials::af::flex_tiny6_double>);
    def("remove_if_not", &remove_if_not <dials::af::flex_miller_index>);
}

} // namespace = boost_python

}} // namespace = dials::array_family
