
#include <boost/python.hpp>
#include <boost/python/def.hpp>

using namespace boost::python;

namespace dials { namespace spot_prediction { 

namespace boost_python {

void export_index_generator();
void export_xds_rotation_angles();
void export_spot_predictor();

BOOST_PYTHON_MODULE(dials_spot_prediction_ext)
{
    export_index_generator();
    export_xds_rotation_angles();
    export_spot_predictor();
}

} // namespace = boost_python

}} // namespace = dials::spot_prediction
