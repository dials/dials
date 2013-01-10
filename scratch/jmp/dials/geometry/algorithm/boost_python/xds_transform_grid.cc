
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../xds_transform_grid.h"

using namespace boost::python;

namespace dials { namespace geometry { namespace algorithm { 
    
namespace boost_python {

void export_xds_transform_grid() 
{
    class_ <xds_transform_grid> ("xds_transform_grid")
        .def(init <std::size_t,
                   scitbx::vec3 <std::size_t>,
                   double,
                   double,
                   double> ((
                arg("n_ref"),
                arg("origin"),
                arg("sigma_divergence"), 
                arg("sigma_mosaicity"),
                arg("n_sigma") = 10)))
        .add_property("n_reflections", 
            &xds_transform_grid::get_n_reflections)
        .add_property("size", 
            &xds_transform_grid::get_size)
        .add_property("origin",
            &xds_transform_grid::get_origin)
        .add_property("step_size", 
            &xds_transform_grid::get_step_size)
        .add_property("sigma_divergence", 
            &xds_transform_grid::get_sigma_divergence)
        .add_property("sigma_mosaicity", 
            &xds_transform_grid::get_sigma_mosaicity)
        .add_property("delta_divergence", 
            &xds_transform_grid::get_delta_divergence)
        .add_property("delta_mosaicity", 
            &xds_transform_grid::get_delta_mosaicity)
        .add_property("data",
            &xds_transform_grid::get_data);       
}

}

}}}
