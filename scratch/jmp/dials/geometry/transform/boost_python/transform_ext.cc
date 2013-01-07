
#include <boost/python.hpp>
#include <boost/python/def.hpp>

using namespace boost::python;

namespace dials { namespace geometry { namespace transform {
    namespace boost_python {

void export_from_hkl_to_beam_vector();
void export_from_beam_vector_to_detector();
void export_from_hkl_to_detector();

BOOST_PYTHON_MODULE(transform_ext)
{
    export_from_hkl_to_beam_vector();
    export_from_beam_vector_to_detector();
    export_from_hkl_to_detector();
}

}}}}
