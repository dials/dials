#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/algorithms/background/lui_2d_background.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_lui_2d_background() {
    def("flat_background_flex_2d", &flat_background_flex_2d, (arg("data2d"), arg("mask2d") ));
    //def("smooth_3d", &smooth_3d, (arg("data2d"),arg("a") = 5) );
  }

}}} // namespace = dials::algorithms::boost_python
