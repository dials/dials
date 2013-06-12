#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/algorithms/background/lui_2d_background.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_lui_find_peak_smoothing() {
    def("bakground_subtract_2d", &bakground_subtract_2d, (arg("data2d") ));
    //def("smooth_3d", &smooth_3d, (arg("data2d"),arg("a") = 5) );
  }

}}} // namespace = dials::algorithms::boost_python
