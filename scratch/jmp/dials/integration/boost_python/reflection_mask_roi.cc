
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../reflection_mask_roi.h"

using namespace boost::python;

namespace dials { namespace integration { 

namespace boost_python {

void export_reflection_mask_roi()
{
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
        .def("calculate",
            &ReflectionMaskRoi::calculate, (
                arg("s1"),
                arg("phi")));
}

} // namespace = boost_python

}} // namespace = dials::integration
