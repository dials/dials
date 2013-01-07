
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../from_hkl_to_beam_vector.h"

using namespace boost::python;
using namespace dials::geometry::transform;

namespace dials { namespace geometry { namespace transform { 
    namespace boost_python {

void export_from_hkl_to_beam_vector() 
{
    class_<from_hkl_to_beam_vector>("from_hkl_to_beam_vector")
        .def(init <scitbx::vec3 <double>, scitbx::vec3 <double>,
                   scitbx::vec3 <double>, scitbx::vec3 <double>,
                   scitbx::vec3 <double> > ((
            arg("b1_star"), arg("b2_star"), arg("b3_star"), 
            arg("s0"), arg("m2"))))
        .def("apply", &from_hkl_to_beam_vector::apply, (
            arg("hkl"), arg("phi")));
}

}}}}
