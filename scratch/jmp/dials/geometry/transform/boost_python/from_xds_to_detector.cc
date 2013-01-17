
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../from_xds_to_detector.h"

using namespace boost::python;

namespace dials { namespace geometry { namespace transform { 
    
namespace boost_python {

void export_from_xds_to_detector() 
{
    class_ <FromXdsToDetector> ("FromXdsToDetector")
        .def(init <FromXdsToBeamVector,
                   FromBeamVectorToDetector> ((
                arg("from_xds_to_beam_vector"), 
                arg("from_beam_vector_to_detector"))))
        .def(init <XdsCoordinateSystem,
                   scitbx::vec3 <double>,
                   DetectorCoordinateSystem,
                   scitbx::vec2 <double>,
                   scitbx::vec2 <double>,
                   double> ((
                arg("xcs"), 
                arg("s1"),
                arg("dcs"),
                arg("pixel_size"),
                arg("origin"),
                arg("distance"))))
        .def("apply", 
            &FromXdsToDetector::apply, (
                arg("c")));
}

} // namespace = boost_python

}}} // namespace = dials::geometry::transform
