
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../from_xds_e3_to_phi.h"

using namespace boost::python;

namespace dials { namespace geometry { namespace transform { 
    
namespace boost_python {

void export_from_xds_e3_to_phi() 
{
    class_ <FromXdsE3ToPhi> ("FromXdsE3ToPhi")
        .def(init <double,
                   double> ((
                arg("zeta"), 
                arg("phi"))))
        .def("apply", 
            &FromXdsE3ToPhi::apply, (
                arg("e3")));
}

} // namespace = boost_python

}}} // namespace = dials::geometry::transform
