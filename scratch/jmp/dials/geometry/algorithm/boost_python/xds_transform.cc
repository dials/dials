
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../xds_transform.h"

using namespace boost::python;

namespace dials { namespace geometry { namespace algorithm { 
    
namespace boost_python {

void export_xds_transform() 
{
    class_ <xds_transform> ("xds_transform")
        .def(init <scitbx::af::flex_int,
                   scitbx::vec3 <int>,
                   scitbx::vec3 <int>,
                   scitbx::vec3 <int>,
                   scitbx::vec3 <double>,
                   double,
                   double,
                   double,
                   double > ((
                arg("image"),
                arg("image_size"),
                arg("roi_size"), 
                arg("grid_size"),
                arg("step_size"),
                arg("starting_frame"),
                arg("starting_angle"),
                arg("oscillation_range"),
                arg("sigma_mosaicity"))))
        .def("calculate", 
            &xds_transform::calculate,(
                arg("frame"),
                arg("phi"),
                arg("zeta")))
        .def("calculate2",
            &xds_transform::calculate2)
        .def("calculate3",
            &xds_transform::calculate3);

    class_ <xds_transform_e3_fraction> ("xds_transform_e3_fraction")
        .def(init <int,
                   int,
                   double,    
                   double,
                   double,
                   double > ((
                arg("roi_size_z"), 
                arg("grid_size_e3"),
                arg("step_size_e3"),
                arg("starting_angle"),
                arg("oscillation_range"),
                arg("sigma_mosaicity"))))
        .def("calculate", 
            &xds_transform_e3_fraction::calculate,(
                arg("frame"),
                arg("phi"),
                arg("zeta")));
}

}

}}}
