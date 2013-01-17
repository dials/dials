
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../from_hkl_to_detector.h"

using namespace boost::python;
using namespace dials::geometry::transform;

namespace dials { namespace geometry { namespace transform { 
    
namespace boost_python {

void export_from_hkl_to_detector() 
{
    class_ <FromHklToDetector> ("FromHklToDetector")
        .def(init <FromHklToBeamVector, 
                   FromBeamVectorToDetector > ((
                arg("hkl_to_s1"), 
                arg("s1_to_xy"))))
        .def(init <ReciprocalLatticeCoordinateSystem,
                   scitbx::vec3 <double>,
                   scitbx::vec3 <double>,
                   DetectorCoordinateSystem,
                   scitbx::vec2 <double>,
                   scitbx::vec2 <double>,
                   double> ((
                arg("rlcs"),
                arg("s0"),
                arg("m2"),
                arg("dcs"),
                arg("pixel_size"),
                arg("origin"),
                arg("distance"))))               
        .def("apply", 
            &FromHklToDetector::apply, (
                arg("hkl"), 
                arg("phi")));
}

} // namespace = boost_python

}}} // namespace = dials::geometry::transform
