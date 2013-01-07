
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../from_beam_vector_to_xds.h"

using namespace boost::python;

namespace dials { namespace geometry { namespace transform { 
    
namespace boost_python {

void export_from_beam_vector_to_xds() 
{
    class_ <from_beam_vector_to_xds> ("from_beam_vector_to_xds")
        .def(init <xds_coordinate_system, 
                   scitbx::vec3 <double>,
                   double> ((
                arg("xcs"), 
                arg("s1"), 
                arg("phi"))))
        .def("apply", 
            &from_beam_vector_to_xds::apply, (
                arg("s_dash"),
                arg("phi_dash")));
}

}

}}}
