
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../outlier_detection/mahalanobis.h"

using namespace boost::python;

namespace dials { namespace refinement { namespace boost_python {

  void export_mahalanobis() {
    def("maha_dist_sq", &maha_dist_sq, (arg("obs"), arg("center"), arg("cov")));
  }

}}}  // namespace dials::refinement::boost_python
