#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/algorithms/peak_finding/lui_find_peak_helper.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_lui_find_peak_helper() {
    def("lui_hi", &hello);

  }

}}} // namespace = dials::algorithms::boost_python
