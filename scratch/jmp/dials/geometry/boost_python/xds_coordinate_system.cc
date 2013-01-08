
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../xds_coordinate_system.h"

using namespace boost::python;

namespace dials { namespace geometry { 

namespace boost_python {

void export_xds_coordinate_system()
{
    class_ <xds_coordinate_system> (
            "xds_coordinate_system")
        .def(init <scitbx::vec3 <double>,
                   scitbx::vec3 <double>,
                   scitbx::vec3 <double>,
                   double > ((
                arg("s0"), 
                arg("s1"), 
                arg("m2"),
                arg("phi"))))
        .add_property("e1",   &xds_coordinate_system::get_e1_axis)
        .add_property("e2",   &xds_coordinate_system::get_e2_axis)
        .add_property("e3",   &xds_coordinate_system::get_e3_axis)
        .add_property("zeta", &xds_coordinate_system::get_zeta);
}

}}}
