
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../from_xds_to_detector.h"

using namespace boost::python;

namespace dials { namespace geometry { namespace transform { 
    
namespace boost_python {

void export_from_xds_to_detector() 
{
    class_ <FromXdsToDetector> ("FromXdsToDetector")
        .def(init <const FromXdsToBeamVector&,
                   const FromBeamVectorToDetector&> ((
                arg("from_xds_to_beam_vector"), 
                arg("from_beam_vector_to_detector"))))
        .def(init <const XdsCoordinateSystem&,
                   scitbx::vec3 <double>,
                   const equipment::Detector&> ((
                arg("xcs"), 
                arg("s1"),
                arg("detector"))))
        .def("apply", 
            &FromXdsToDetector::apply, (
                arg("c")));
}

} // namespace = boost_python

}}} // namespace = dials::geometry::transform
