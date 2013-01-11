
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../from_beam_vector_to_detector.h"

using namespace boost::python;
using namespace dials::geometry::transform;

namespace dials { namespace geometry { namespace transform { 
    
namespace boost_python {

void export_from_beam_vector_to_detector() 
{
    class_ <FromBeamVectorToDetector> ("FromBeamVectorToDetector")
        .def(init <DetectorCoordinateSystem, 
                   scitbx::vec2 <double>,
                   double> ((
                arg("dcs"), 
                arg("origin"), 
                arg("distance"))))
        .def("apply", 
            &FromBeamVectorToDetector::apply, (
                arg("s1")));
}

} // namespace = boost_python

}}} // namespace = dials::geometry::transform
