
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../angle_filter.h"

using namespace boost::python;

namespace dials { namespace spot_prediction { 

namespace boost_python {

void export_angle_filter()
{
    def("angle_filter", &angle_filter, (
        arg("angle"),
        arg("range"),
        arg("deg") = false));
}

} // namespace = boost_python

}} // namespace = dials::spot_prediction
