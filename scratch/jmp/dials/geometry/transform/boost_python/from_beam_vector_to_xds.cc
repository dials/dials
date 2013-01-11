
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../from_beam_vector_to_xds.h"

using namespace boost::python;

namespace dials { namespace geometry { namespace transform { 
    
namespace boost_python {

void export_from_beam_vector_to_xds() 
{
    class_ <FromBeamVectorToXds> ("FromBeamVectorToXds")
        .def(init <XdsCoordinateSystem, 
                   scitbx::vec3 <double>,
                   double> ((
                arg("xcs"), 
                arg("s1"), 
                arg("phi"))))
        .def("apply", 
            &FromBeamVectorToXds::apply, (
                arg("s_dash"),
                arg("phi_dash")));
}

} // namespace = boost_python

}}} // namespace = dials::geometry::transform
