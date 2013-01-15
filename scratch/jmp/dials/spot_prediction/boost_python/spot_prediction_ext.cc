
#include <boost/python.hpp>
#include <boost/python/def.hpp>

using namespace boost::python;

namespace dials { namespace spot_prediction { 

namespace boost_python {

void export_reflection_mask();
void export_background_intensity();
void export_subtract_background();

BOOST_PYTHON_MODULE(dials_spot_prediction_ext)
{
    export_reflection_mask();
    export_background_intensity();
    export_subtract_background();
}

} // namespace = boost_python

}} // namespace = dials::spot_prediction
