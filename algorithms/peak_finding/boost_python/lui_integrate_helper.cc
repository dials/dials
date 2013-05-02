#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/algorithms/peak_finding/lui_integrate_helper.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_lui_integrate_helper() {
//    def("lui_hi", &hello);
    def("ref_2d", &ref_2d, (arg("data2d"),arg("a") = 5) );
  }

}}} // namespace = dials::algorithms::boost_python
