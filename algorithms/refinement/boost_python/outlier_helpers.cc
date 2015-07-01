
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../outlier_helpers.h"

using namespace boost::python;

namespace dials { namespace refinement { namespace boost_python {

  void export_outlier_helpers()
  {
    def("mcd_consistency", &mcd_consistency, (
      arg("p"),
      arg("alpha")));
  }

}}} // namespace dials::refinement::boost_python
