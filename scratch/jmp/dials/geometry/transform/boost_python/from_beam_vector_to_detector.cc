
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../from_beam_vector_to_detector.h"

using namespace boost::python;
using namespace dials::geometry::transform;

namespace dials { namespace geometry { namespace transform { 
    
namespace boost_python {

void export_from_beam_vector_to_detector() 
{
    // Apply to a single beam vector
    scitbx::vec2 <double> (FromBeamVectorToDetector::*apply_single)(
        scitbx::vec3 <double>) const = &FromBeamVectorToDetector::apply;

    // Apply to array of beam vectors
    flex_vec2_double (FromBeamVectorToDetector::*apply_array)(
        const flex_vec3_double&, scitbx::af::flex_bool&) const = 
            &FromBeamVectorToDetector::apply;

    class_ <FromBeamVectorToDetector> ("FromBeamVectorToDetector")
        .def(init <const equipment::Detector&> ((
                arg("detector"))))
        .def("apply", apply_single, (
                arg("s1")))
        .def("apply", apply_array, (
                arg("s1"),
                arg("status")));
}

} // namespace = boost_python

}}} // namespace = dials::geometry::transform
