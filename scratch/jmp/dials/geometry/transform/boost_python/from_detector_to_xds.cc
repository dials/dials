
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../from_detector_to_xds.h"

using namespace boost::python;

namespace dials { namespace geometry { namespace transform { 
    
namespace boost_python {

void export_from_detector_to_xds() 
{
    class_ <from_detector_to_xds> ("from_detector_to_xds")
        .def(init <from_detector_to_beam_vector, 
                   from_beam_vector_to_xds,
                   double> ((
                arg("xy_to_s1"), 
                arg("s1_to_xds"),
                arg("wavelength"))))
        .def(init <detector_coordinate_system,
                   scitbx::vec2 <double>,
                   double,
                   xds_coordinate_system,
                   scitbx::vec3 <double>,
                   double,
                   double> ((
                arg("dcs"),
                arg("origin"),
                arg("distance"),
                arg("s1"),
                arg("phi"),
                arg("wavelength"))))        
        .def("apply", 
            &from_detector_to_xds::apply, (
                arg("xy"),
                arg("phi_dash")));
}

}

}}}
