
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <scitbx/array_family/flex_types.h>
#include "../centroid.h"

using namespace boost::python;

namespace dials { namespace integration { 

namespace boost_python {

void export_centroid()
{
    scitbx::vec2 <double> (*centroid2d_int)(const scitbx::af::flex_int&, 
        scitbx::af::tiny <int, 4>) = &centroid2d <scitbx::af::flex_int>;

    scitbx::vec2 <double> (*centroid2d_int_mask)(const scitbx::af::flex_int&, 
        const scitbx::af::flex_int&, scitbx::af::tiny <int, 4>, int) = 
            &centroid2d <scitbx::af::flex_int, scitbx::af::flex_int>;

    scitbx::vec3 <double> (*centroid3d_int)(const scitbx::af::flex_int&, 
        scitbx::af::tiny <int, 6>) = &centroid3d <scitbx::af::flex_int>;

    scitbx::vec3 <double> (*centroid3d_int_mask)(const scitbx::af::flex_int&, 
        const scitbx::af::flex_int&, scitbx::af::tiny <int, 6>, int) = 
            &centroid3d <scitbx::af::flex_int, scitbx::af::flex_int>;

    scitbx::vec3 <double> (*centroid_reflection_int)(const scitbx::af::flex_int&,
        const scitbx::af::flex_int&, scitbx::af::tiny <int, 6>, int, double) =
            &centroid_reflection <scitbx::af::flex_int, scitbx::af::flex_int>;

    def("centroid2d", centroid2d_int);
    def("centroid3d", centroid3d_int);
    def("centroid2d", centroid2d_int_mask);
    def("centroid3d", centroid3d_int_mask);
    def("centroid_reflection", centroid_reflection_int);
}

} // namespace = boost_python

}} // namespace = dials::spot_prediction
