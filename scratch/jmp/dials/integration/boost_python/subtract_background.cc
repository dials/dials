
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../subtract_background.h"

using namespace boost::python;

namespace dials { namespace integration { 

namespace boost_python {

void export_subtract_background()
{
    class_ <SubtractBackground> ("SubtractBackground")
        .def(init <scitbx::af::flex_int,
                   scitbx::vec3 <int>,
                   double,
                   double> ((
            arg("image_volume"),
            arg("roi_size"),
            arg("delta"),
            arg("max_iter"))))
        .def("subtract", &SubtractBackground::subtract, (
            arg("xyz")));
}

} // namespace = boost_python

}} // namespace = dials::spot_prediction
