
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <cctbx/miller.h>
#include "../coordinate_filter.h"

using namespace boost::python;

namespace dials { namespace spot_prediction { 

namespace boost_python {

typedef cctbx::miller::index <> miller_index;
typedef scitbx::af::flex <scitbx::vec2 <int> >::type flex_vec2_int;
typedef scitbx::af::flex <scitbx::vec3 <int> >::type flex_vec3_int;
typedef scitbx::af::flex <miller_index >::type flex_miller_index;

void export_coordinate_filter()
{
    def("in_range", &in_range, (
        arg("x"),
        arg("range")));
    def("in_rect", &in_rect, (
        arg("xy"),
        arg("xrange"),
        arg("yrange")));
    def("in_volume", &in_volume, (
        arg("xyz"),
        arg("xrange"),
        arg("yrange"),
        arg("zrange")));

    def("remove_if", &remove_if <scitbx::af::flex_bool>);
    def("remove_if", &remove_if <scitbx::af::flex_int>);
    def("remove_if", &remove_if <scitbx::af::flex_long>);
    def("remove_if", &remove_if <scitbx::af::flex_size_t>);
    def("remove_if", &remove_if <scitbx::af::flex_float>);
    def("remove_if", &remove_if <scitbx::af::flex_double>);
    def("remove_if", &remove_if <scitbx::af::flex_complex_double>);
    def("remove_if", &remove_if <flex_vec2_int>);
    def("remove_if", &remove_if <flex_vec2_double>);
    def("remove_if", &remove_if <flex_vec3_int>);
    def("remove_if", &remove_if <flex_vec3_double>);
    def("remove_if", &remove_if <flex_miller_index>);

    def("remove_if_not", &remove_if_not <scitbx::af::flex_bool>);
    def("remove_if_not", &remove_if_not <scitbx::af::flex_int>);
    def("remove_if_not", &remove_if_not <scitbx::af::flex_long>);
    def("remove_if_not", &remove_if_not <scitbx::af::flex_size_t>);
    def("remove_if_not", &remove_if_not <scitbx::af::flex_float>);
    def("remove_if_not", &remove_if_not <scitbx::af::flex_double>);
    def("remove_if_not", &remove_if_not <scitbx::af::flex_complex_double>);
    def("remove_if_not", &remove_if_not <flex_vec2_int>);
    def("remove_if_not", &remove_if_not <flex_vec2_double>);
    def("remove_if_not", &remove_if_not <flex_vec3_int>);
    def("remove_if_not", &remove_if_not <flex_vec3_double>);
    def("remove_if_not", &remove_if_not <flex_miller_index>);
}

} // namespace = boost_python

}} // namespace = dials::spot_prediction
