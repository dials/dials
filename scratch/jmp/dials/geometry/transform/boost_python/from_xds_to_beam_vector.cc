
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../from_xds_to_beam_vector.h"

using namespace boost::python;

namespace dials { namespace geometry { namespace transform { 
    
namespace boost_python {

void export_from_xds_to_beam_vector() 
{
    class_ <FromXdsToBeamVector> ("FromXdsToBeamVector")
        .def(init <XdsCoordinateSystem, 
                   scitbx::vec3 <double> > ((
                arg("xcs"), 
                arg("s1"))))
        .def("apply", 
            &FromXdsToBeamVector::apply, (
                arg("c")));
}

} // namespace = boost_python

}}} // namespace = dials::geometry::transform
