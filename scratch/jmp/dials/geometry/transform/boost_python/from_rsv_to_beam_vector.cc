
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../from_rsv_to_beam_vector.h"

using namespace boost::python;
using namespace dials::geometry::transform;

namespace dials { namespace geometry { namespace transform { 
    
namespace boost_python {

void export_from_rsv_to_beam_vector() 
{
    // Apply to a single beam vector
    scitbx::vec3 <double> (FromRsvToBeamVector::*apply_single)(
        scitbx::vec3 <double>, double) = &FromRsvToBeamVector::apply;

    // Apply to array of beam vectors
    flex_vec3_double (FromRsvToBeamVector::*apply_array)(
        const flex_vec3_double&, const scitbx::af::flex_double&) = 
            &FromRsvToBeamVector::apply;

    class_ <FromRsvToBeamVector> ("FromRsvToBeamVector")
        .def(init <const equipment::Beam&,
                   const equipment::Goniometer&> ((
                arg("beam"),
                arg("goniometer"))))
        .def("apply", apply_single, (
                arg("reciprocal_space_vector"), 
                arg("rotation_angle")))
        .def("apply", apply_array, (
                arg("reciprocal_space_vectors"),
                arg("rotation_angles")));
}

} // namespace = boost_python

}}} // namespace = dials::geometry::transform
