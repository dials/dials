
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../xds_rotation_angles.h"

using namespace boost::python;

namespace dials { namespace spot_prediction { 

namespace boost_python {

void export_xds_rotation_angles()
{
    scitbx::vec2 <double> (XdsRotationAngles::*calculate_pstar0) (
        scitbx::vec3 <double>) const = &XdsRotationAngles::calculate;

    scitbx::vec2 <double> (XdsRotationAngles::*calculate_miller) (
        cctbx::miller::index <>, scitbx::mat3 <double>) const = 
            &XdsRotationAngles::calculate;

    class_ <XdsRotationAngles> ("XdsRotationAngles", no_init)
        .def(init <scitbx::vec3 <double>, 
                   scitbx::vec3 <double> > ((
            arg("beam_direction"),
            arg("rotation_axis"))))
        .def("calculate", calculate_pstar0)
        .def("calculate", calculate_miller);
}

} // namespace = boost_python

}} // namespace = dials::spot_prediction
