
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../outlier_detection/outlier_helpers.h"

using namespace boost::python;

namespace dials { namespace refinement { namespace boost_python {

  void export_outlier_helpers() {
    def("qchisq", &qchisq, (arg("p"), arg("df")));
    def("mcd_consistency", &mcd_consistency, (arg("df"), arg("alpha")));
  }

}}}  // namespace dials::refinement::boost_python
