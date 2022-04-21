
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../gallego_yezzi.h"

using namespace boost::python;

namespace dials { namespace refinement { namespace boost_python {

  void export_gallego_yezzi() {
    def("dRq_de", &dRq_de, (arg("theta"), arg("e1"), arg("q")));
  }

}}}  // namespace dials::refinement::boost_python
