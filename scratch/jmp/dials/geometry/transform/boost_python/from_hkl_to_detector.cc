
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
        .def(init <const FromHklToBeamVector&, 
                   const FromBeamVectorToDetector&> ((
                arg("hkl_to_s1"), 
                arg("s1_to_xy"))))
        .def(init <scitbx::mat3 <double>,
                   scitbx::vec3 <double>,
                   scitbx::vec3 <double>,
                   const equipment::Detector&> ((
                arg("ub_matrix"),
                arg("s0"),
                arg("m2"),
                arg("detector"))))               
        .def("apply", 
            &FromHklToDetector::apply, (
                arg("hkl"), 
                arg("phi")));
}

} // namespace = boost_python

}}} // namespace = dials::geometry::transform
