
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../reflection_mask_roi.h"

using namespace boost::python;

namespace dials { namespace integration { 

namespace boost_python {

void export_reflection_mask_roi()
{
    scitbx::af::tiny <int, 6> (ReflectionMaskRoi::*calculate_single)(
        scitbx::vec3 <double>, double) = &ReflectionMaskRoi::calculate;

    af::flex_tiny6_int (ReflectionMaskRoi::*calculate_array) (
        const af::flex_vec3_double &, const scitbx::af::flex_double &) =
            &ReflectionMaskRoi::calculate;

    class_ <ReflectionMaskRoi> ("ReflectionMaskRoi")
        .def(init <const equipment::Beam &,
                   const equipment::Detector &,
                   const equipment::Goniometer &,
                   double,
                   double > ((
            arg("beam"), 
            arg("detector"),
            arg("goniometer"),
            arg("delta_divergence"),
            arg("delta_mosaicity"))))
        .def("calculate", calculate_single, (
                arg("s1"),
                arg("phi")))
        .def("calculate", calculate_array, (
                arg("s1"),
                arg("phi")));
}

} // namespace = boost_python

}} // namespace = dials::integration
