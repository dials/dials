
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../rotation_angles.h"

using namespace boost::python;

namespace dials { namespace spot_prediction { 

namespace boost_python {

void export_rotation_angles()
{
    class_ <RotationAngles> ("RotationAngles", no_init)
        .def(init <double, 
                   scitbx::mat3 <double>,
                   double, 
                   scitbx::vec3 <double> > ((
            arg("d_min"),
            arg("ub_matrix"),
            arg("wavelength"),
            arg("rotation_axis"))))
//            arg("rotation_angle_range"))))
        .def("calculate", &RotationAngles::calculate)
        .def("array_indices", &RotationAngles::get_array_indices)
        .def("miller_indices", &RotationAngles::get_miller_indices);
}

} // namespace = boost_python

}} // namespace = dials::spot_prediction
