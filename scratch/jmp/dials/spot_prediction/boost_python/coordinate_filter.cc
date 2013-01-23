
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../coordinate_filter.h"

using namespace boost::python;

namespace dials { namespace spot_prediction { 

namespace boost_python {

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
}

} // namespace = boost_python

}} // namespace = dials::spot_prediction
