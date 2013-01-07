
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../from_detector_to_beam_vector.h"

using namespace boost::python;
using namespace dials::geometry::transform;

namespace dials { namespace geometry { namespace transform { 
    
namespace boost_python {

void export_from_detector_to_beam_vector() 
{
    class_ <from_detector_to_beam_vector> ("from_detector_to_beam_vector")
        .def(init <detector_coordinate_system, 
                   scitbx::vec2 <double>,
                   double> ((
                arg("dcs"), 
                arg("origin"), 
                arg("distance"))))
        .def("apply", 
            &from_detector_to_beam_vector::apply, (
                arg("xy")));
}

}

}}}
