
#include <boost/python.hpp>
#include <boost/python/def.hpp>

using namespace boost::python;

namespace dials { namespace geometry { namespace transform {

namespace boost_python {

void export_from_hkl_to_beam_vector();
void export_from_hkl_to_detector();
void export_from_beam_vector_to_detector();
void export_from_beam_vector_to_xds();
void export_from_detector_to_beam_vector();
void export_from_detector_to_xds();
void export_from_xds_to_beam_vector();
void export_from_xds_to_detector();
void export_from_xds_e3_to_phi();
void export_from_beam_vector_to_image_volume();

BOOST_PYTHON_MODULE(dials_geometry_transform_ext)
{
    export_from_hkl_to_beam_vector();
    export_from_hkl_to_detector();
    export_from_beam_vector_to_detector();
    export_from_beam_vector_to_xds();
    export_from_detector_to_beam_vector();
    export_from_detector_to_xds();
    export_from_xds_to_beam_vector();
    export_from_xds_to_detector();
    export_from_xds_e3_to_phi();
    export_from_beam_vector_to_image_volume();
}

} // namespace = boost_python

}}} // namespace = dials::geometry::transform
