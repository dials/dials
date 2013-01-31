
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../subtract_background.h"

using namespace boost::python;

namespace dials { namespace integration { 

namespace boost_python {

void export_subtract_background()
{
    double (SubtractBackground::*subtract_single)(int, scitbx::af::tiny <int, 6>) =
        &SubtractBackground::subtract;

    scitbx::af::flex_bool (SubtractBackground::*subtract_array)(
        ReflectionList&) = &SubtractBackground::subtract;

    class_ <SubtractBackground> ("SubtractBackground")
        .def(init <scitbx::af::flex_int,
                   scitbx::af::flex_int,
                   double,
                   double,
                   int> ((
            arg("image_volume"),
            arg("reflection_mask"),
            arg("delta") = 0.1,
            arg("max_iter") = 0.1,
            arg("min_pixels") = 10)))
        .def("subtract", subtract_single, (
            arg("index"),
            arg("roi")))
        .def("subtract", subtract_array, (
            arg("reflections")))
        .def("set_non_reflection_value", 
            &SubtractBackground::set_non_reflection_value);
}

} // namespace = boost_python

}} // namespace = dials::spot_prediction
