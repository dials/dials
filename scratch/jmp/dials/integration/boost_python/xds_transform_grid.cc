
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../xds_transform_grid.h"

using namespace boost::python;

namespace dials { namespace integration { 
    
namespace boost_python {

void export_xds_transform_grid() 
{
    class_ <XdsTransformGrid> ("XdsTransformGrid")
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
            &XdsTransformGrid::get_n_reflections)
        .add_property("size", 
            &XdsTransformGrid::get_size)
        .add_property("origin",
            &XdsTransformGrid::get_origin)
        .add_property("step_size", 
            &XdsTransformGrid::get_step_size)
        .add_property("sigma_divergence", 
            &XdsTransformGrid::get_sigma_divergence)
        .add_property("sigma_mosaicity", 
            &XdsTransformGrid::get_sigma_mosaicity)
        .add_property("delta_divergence", 
            &XdsTransformGrid::get_delta_divergence)
        .add_property("delta_mosaicity", 
            &XdsTransformGrid::get_delta_mosaicity)
        .add_property("data",
            &XdsTransformGrid::get_data);       
}

} // namespace = boost_python

}} // namespace = dials::integration
