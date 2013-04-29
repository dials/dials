#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/algorithms/peak_finding/lui_find_peak_smoothing.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_lui_find_peak_smoothing() {
    def("smooth_2d", &smooth_2d, (arg("data2d"),arg("a") = 5) );
    def("smooth_3d", &smooth_3d, (arg("data2d"),arg("a") = 5) );
  }

}}} // namespace = dials::algorithms::boost_python
