
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../from_detector_to_xds.h"

using namespace boost::python;

namespace dials { namespace geometry { namespace transform { 
    
namespace boost_python {

void export_from_detector_to_xds() 
{
    class_ <FromDetectorToXds> ("FromDetectorToXds")
        .def(init <FromDetectorToBeamVector, 
                   FromBeamVectorToXds,
                   double> ((
                arg("xy_to_s1"), 
                arg("s1_to_xds"),
                arg("wavelength"))))
        .def(init <DetectorCoordinateSystem,
                   scitbx::vec2 <double>,
                   scitbx::vec2 <double>,
                   double,
                   XdsCoordinateSystem,
                   scitbx::vec3 <double>,
                   double,
                   double> ((
                arg("dcs"),
                arg("pixel_size"),
                arg("origin"),
                arg("distance"),
                arg("s1"),
                arg("phi"),
                arg("wavelength"))))        
        .def("apply", 
            &FromDetectorToXds::apply, (
                arg("xy"),
                arg("phi_dash")));
}

} // namespace = boost_python

}}} // namespace = dials::geometry::transform
