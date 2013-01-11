
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../xds_coordinate_system.h"

using namespace boost::python;

namespace dials { namespace geometry { 

namespace boost_python {

void export_xds_coordinate_system()
{
    class_ <XdsCoordinateSystem> ("XdsCoordinateSystem")
        .def(init <scitbx::vec3 <double>,
                   scitbx::vec3 <double>,
                   scitbx::vec3 <double>,
                   double > ((
                arg("s0"), 
                arg("s1"), 
                arg("m2"),
                arg("phi"))))
        .add_property("e1",   &XdsCoordinateSystem::get_e1_axis)
        .add_property("e2",   &XdsCoordinateSystem::get_e2_axis)
        .add_property("e3",   &XdsCoordinateSystem::get_e3_axis)
        .add_property("zeta", &XdsCoordinateSystem::get_zeta);
}

} // namespace = boost_python

}} // namespace = dials::geometry
