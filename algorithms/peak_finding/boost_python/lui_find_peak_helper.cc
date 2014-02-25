#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/algorithms/peak_finding/lui_find_peak_helper.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_lui_find_peak_helper() {
    def("find_mask_2d", &find_mask_2d, (
    arg("data2d")
    ,arg("data2dsmoth")
    ,arg("a") = 5

    )

    );
  }

}}}
