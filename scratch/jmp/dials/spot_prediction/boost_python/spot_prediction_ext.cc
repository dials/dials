
#include <boost/python.hpp>
#include <boost/python/def.hpp>

using namespace boost::python;

namespace dials { namespace spot_prediction { 

namespace boost_python {

void export_reflection_mask();
void export_background_intensity();
void export_subtract_background();
void export_index_generator();
void export_rotation_angles();
void export_angle_filter();
void export_coordinate_filter();

BOOST_PYTHON_MODULE(dials_spot_prediction_ext)
{
    export_reflection_mask();
    export_background_intensity();
    export_subtract_background();
    export_index_generator();
    export_rotation_angles();
    export_angle_filter();
    export_coordinate_filter();
}

} // namespace = boost_python

}} // namespace = dials::spot_prediction
