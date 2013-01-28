
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <scitbx/array_family/flex_types.h>
#include "../centroid.h"

using namespace boost::python;

namespace dials { namespace integration { 

namespace boost_python {

void export_centroid()
{
    def("centroid2d", &centroid2d <scitbx::af::flex_int>);
    def("centroid3d", &centroid3d <scitbx::af::flex_int>);
    def("centroid2d", &centroid2d <scitbx::af::flex_int, scitbx::af::flex_int>);
    def("centroid3d", &centroid3d <scitbx::af::flex_int, scitbx::af::flex_int>);
}

} // namespace = boost_python

}} // namespace = dials::spot_prediction
