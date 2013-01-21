
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../from_hkl_to_rsv.h"

using namespace boost::python;
using namespace dials::geometry::transform;

namespace dials { namespace geometry { namespace transform { 
    
namespace boost_python {

void export_from_hkl_to_rsv() 
{
    // Apply to a single beam vector
    scitbx::vec3 <double> (FromHklToRsv::*apply_single)(
        miller_index, double) = &FromHklToRsv::apply;

    // Apply to array of beam vectors
    flex_vec3_double (FromHklToRsv::*apply_array)(
        flex_miller_index, scitbx::af::flex_double) = &FromHklToRsv::apply;

    class_ <FromHklToRsv> ("FromHklToRsv")
        .def(init <scitbx::mat3 <double>,
                   scitbx::vec3 <double> > ((
                arg("ub"),
                arg("m2"))))
        .def("apply", apply_single, (
                arg("h"),
                arg("phi")))
        .def("apply", apply_array, (
                arg("h"),
                arg("phi")));
}

} // namespace = boost_python

}}} // namespace = dials::geometry::transform
