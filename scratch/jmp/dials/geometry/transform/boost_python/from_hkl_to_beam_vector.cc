
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../from_hkl_to_beam_vector.h"

using namespace boost::python;

namespace dials { namespace geometry { namespace transform { 
    
namespace boost_python {

void export_from_hkl_to_beam_vector() 
{
    // Apply to a sinfle hkl, phi pair
    scitbx::vec3 <double> (FromHklToBeamVector::*apply_single)(
        miller_index, double) = &FromHklToBeamVector::apply;

    // Apply to array of hkl, phi pairs
    flex_vec3_double (FromHklToBeamVector::*apply_array)(
        flex_miller_index, scitbx::af::flex_double) = &FromHklToBeamVector::apply;

    class_ <FromHklToBeamVector> ("FromHklToBeamVector")
        .def(init <scitbx::mat3 <double>, 
                   scitbx::vec3 <double>,
                   scitbx::vec3 <double> > ((
                arg("ub_matrix"), 
                arg("s0"), 
                arg("m2"))))
        .def("apply", apply_single, (
                arg("hkl"), 
                arg("phi")))
        .def("apply", apply_array, (
                arg("hkl"), 
                arg("phi")));
}

} // namespace = boost_python

}}} // namespace = dials::geometry::transform
