
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../from_hkl_to_beam_vector.h"

using namespace boost::python;

namespace dials { namespace geometry { namespace transform { 
    
namespace boost_python {

void export_from_hkl_to_beam_vector() 
{
    class_ <FromHklToBeamVector> ("FromHklToBeamVector")
        .def(init <ReciprocalLatticeCoordinateSystem, 
                   scitbx::vec3 <double>,
                   scitbx::vec3 <double> > ((
                arg("rlcs"), 
                arg("s0"), 
                arg("m2"))))
        .def("apply", &FromHklToBeamVector::apply, (
                arg("hkl"), 
                arg("phi")));
}

} // namespace = boost_python

}}} // namespace = dials::geometry::transform
