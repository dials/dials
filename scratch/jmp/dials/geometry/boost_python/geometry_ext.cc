
#include <boost/python.hpp>
#include <boost/python/def.hpp>

using namespace boost::python;

namespace dials { namespace geometry { namespace boost_python {

void export_detector_coordinate_system();
void export_reciprocal_lattice_coordinate_system();

BOOST_PYTHON_MODULE(geometry_ext)
{
    export_detector_coordinate_system();
    export_reciprocal_lattice_coordinate_system();
}

}}}
