
#include <boost/python.hpp>
#include <boost/python/def.hpp>

using namespace boost::python;

namespace dials { namespace equipment { 

namespace boost_python {

void export_beam();
void export_detector();
void export_goniometer();

BOOST_PYTHON_MODULE(equipment_ext)
{
    export_beam();
    export_detector();
    export_goniometer();
}

}}

}
