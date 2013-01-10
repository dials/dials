
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../xds_transform.h"

using namespace boost::python;

namespace dials { namespace geometry { namespace algorithm { 
    
namespace boost_python {

void export_xds_transform() 
{
    class_ <xds_transform> ("xds_transform")
        .def(init <xds_transform_grid,
                   scitbx::af::flex_int,
                   scitbx::vec3 <int>,
                   //scitbx::vec3 <int>,
                   //scitbx::vec3 <double>,
                   //double,
                   //double,
                   //double,
                   //double,
                   equipment::detector,
                   equipment::beam,
                   equipment::goniometer,
                   scitbx::vec3 <int>,
                   int > ((
                   //scitbx::vec2 <double>,
                   //detector_coordinate_system,
                   //double,
                   //double,
                   //scitbx::vec3 <double>,
                   //scitbx::vec3 <double> > ((
                arg("xds_transform_grid"),
                arg("image"),
                arg("image_size"),
                //arg("grid_size"),
                //arg("step_size"),
                //arg("starting_frame"),
                //arg("starting_angle"),
                //arg("oscillation_range"),
                //arg("sigma_mosaicity"),
                arg("detector"),
                arg("beam"),
                arg("goniometer"),
                arg("roi_size") = scitbx::vec3 <int> (4, 4, 1), 
                arg("n_div") = 5)))
                //arg("origin"),
                //arg("detector_coordinate_system"),
                //arg("distance"),
                //arg("wavelength"),
                //arg("incident_beam_vector"),
                //arg("rotation_axis"))))    
        .def("calculate",
            &xds_transform::calculate);

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
