
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../reciprocal_lattice_coordinate_system.h"

using namespace boost::python;

namespace dials { namespace geometry { 

namespace boost_python {

void export_reciprocal_lattice_coordinate_system()
{
    class_ <reciprocal_lattice_coordinate_system> (
            "reciprocal_lattice_coordinate_system")
        .def(init <scitbx::vec3 <double>,
                   scitbx::vec3 <double>,
                   scitbx::vec3 <double> > ((
                arg("b1_star"), 
                arg("b2_star"), 
                arg("b3_star"))))
        .def(init <scitbx::mat3 <double> > ((
                arg("ub"))))
        .add_property("b1_star",  
            &reciprocal_lattice_coordinate_system::get_b1_star_axis,
            &reciprocal_lattice_coordinate_system::set_b1_star_axis)
        .add_property("b2_star", 
            &reciprocal_lattice_coordinate_system::get_b2_star_axis,
            &reciprocal_lattice_coordinate_system::set_b2_star_axis)
        .add_property("b3_star", 
            &reciprocal_lattice_coordinate_system::get_b3_star_axis,
            &reciprocal_lattice_coordinate_system::set_b3_star_axis)
        .def("from_ub_matrix",
            &reciprocal_lattice_coordinate_system::from_ub_matrix, (
                arg("ub")))
        .def("to_ub_matrix",
            &reciprocal_lattice_coordinate_system::to_ub_matrix);
 
 }

}}}
