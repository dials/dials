
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../background_intensity.h"

using namespace boost::python;

namespace dials { namespace integration { 

namespace boost_python {

void export_background_intensity()
{
    class_ <BackgroundIntensity> ("BackgroundIntensity")
        .def(init <double, double> ((
            arg("delta") = 0.1,
            arg("max_iter_frac") = 0.1)))
        .def("calculate",
            &BackgroundIntensity::calculate, (
            arg("data")));
}

} // namespace = boost_python

}} // namespace = dials::spot_prediction
