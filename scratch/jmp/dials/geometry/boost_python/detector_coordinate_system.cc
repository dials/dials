
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../detector_coordinate_system.h"

using namespace boost::python;

namespace dials { namespace geometry { 

namespace boost_python {

void export_detector_coordinate_system()
{
    class_ <DetectorCoordinateSystem> ("DetectorCoordinateSystem")
        .def(init <scitbx::vec3 <double>,
                   scitbx::vec3 <double> > ((
                arg("x_axis"), 
                arg("y_axis"))))
        .def(init <scitbx::vec3 <double>,
                   scitbx::vec3 <double>,
                   scitbx::vec3 <double> > ((
                arg("x_axis"), 
                arg("y_axis"), 
                arg("normal"))))
        .add_property("x_axis",  
            &DetectorCoordinateSystem::get_x_axis,
            &DetectorCoordinateSystem::set_x_axis)
        .add_property("y_axis", 
            &DetectorCoordinateSystem::get_y_axis,
            &DetectorCoordinateSystem::set_y_axis)
        .add_property("normal", 
            &DetectorCoordinateSystem::get_normal,
            &DetectorCoordinateSystem::set_normal)
        .def("in_pixel_units",
            &DetectorCoordinateSystem::in_pixel_units, (
                arg("pixel_size")))
        .def("in_si_units",
            &DetectorCoordinateSystem::in_si_units, (
                arg("pixel_size")));     
}

} // namespace = boost_python

}} // namespace = dials::geometry
