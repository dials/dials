#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/algorithms/background/lui_2d_background.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_lui_2d_background() {
    def("flat_background_flex_2d", &flat_background_flex_2d,
        (arg("data2d"), arg("mask2d") ));
    def("curved_background_flex_2d", &curved_background_flex_2d,
        (arg("data2d"), arg("mask2d") ));
    def("get_plane_background_syml_sys_2d", &get_plane_background_syml_sys_2d,
        (arg("data2d"), arg("mask2d"), arg("mat_a"), arg("vec_b") ));
    def("inclined_plane_background_flex_2d", &inclined_plane_background_flex_2d,
        (arg("data2d"), arg("mask2d"), arg("abc_plane"), arg("background2d")));
  }

}}} // namespace = dials::algorithms::boost_python
