
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/algorithms/integration/from_beam_vector_to_xds.h>

using namespace boost::python;

namespace dials { namespace algorithms { namespace boost_python {

  void export_from_beam_vector_to_xds() 
  {
    class_ <FromBeamVectorToXds> ("FromBeamVectorToXds", no_init)
      .def(init <XdsCoordinateSystem, 
                 vec3 <double>,
                 double> ((
          arg("xcs"), 
          arg("s1"), 
          arg("phi"))))
      .def("apply", 
          &FromBeamVectorToXds::apply, (
          arg("s_dash"),
          arg("phi_dash")));
  }

}}} // namespace = dials::algorithms::boost_python
