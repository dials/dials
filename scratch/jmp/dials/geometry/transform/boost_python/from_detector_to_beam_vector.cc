
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../from_detector_to_beam_vector.h"

using namespace boost::python;
using namespace dials::geometry::transform;

namespace dials { namespace geometry { namespace transform { 
    
namespace boost_python {

void export_from_detector_to_beam_vector() 
{
    class_ <FromDetectorToBeamVector> ("FromDetectorToBeamVector")
        .def(init <equipment::Detector> ((
                arg("detector"))))
        .def("apply", 
            &FromDetectorToBeamVector::apply, (
                arg("xy")));
}

} // namespace = boost_python

}}} // namespace = dials::geometry::transform
