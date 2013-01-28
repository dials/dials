
#include <boost/python.hpp>
#include <boost/python/def.hpp>

using namespace boost::python;

namespace dials { namespace integration { 

namespace boost_python {

void export_reflection_mask_roi();
void export_reflection_mask();
void export_background_intensity();
void export_subtract_background();
void export_xds_transform();
void export_xds_transform_grid();
void export_centroid();

BOOST_PYTHON_MODULE(dials_integration_ext)
{
    export_reflection_mask_roi();
    export_reflection_mask();
    export_background_intensity();
    export_subtract_background();
    export_xds_transform();
    export_xds_transform_grid();    
    export_centroid();
}

} // namespace = boost_python

}} // namespace = dials::integration
