
#include <boost/python.hpp>
#include <boost/python/def.hpp>

using namespace boost::python;

namespace dials { namespace geometry { namespace algorithm {

namespace boost_python {

void export_xds_transform();

BOOST_PYTHON_MODULE(algorithm_ext)
{
    export_xds_transform();
}

}

}}}
