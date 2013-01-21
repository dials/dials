
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../angle_filter.h"

using namespace boost::python;

namespace dials { namespace spot_prediction { 

namespace boost_python {

void export_angle_filter()
{
    bool (*is_angle_in_range_single)(double, scitbx::vec2 <double>, bool) =
        &is_angle_in_range;    

    scitbx::af::flex_bool (*is_angle_in_range_array)(
        scitbx::af::flex_double, scitbx::vec2 <double>, bool) =
            &is_angle_in_range;

    def("is_angle_in_range", is_angle_in_range_single, (
        arg("angle"),
        arg("range"),
        arg("deg") = false));

    def("is_angle_in_range", is_angle_in_range_array, (
        arg("angle"),
        arg("range"),
        arg("deg") = false));
}

} // namespace = boost_python

}} // namespace = dials::spot_prediction
