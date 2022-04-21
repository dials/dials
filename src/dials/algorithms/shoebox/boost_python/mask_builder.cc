#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/algorithms/shoebox/mask_builder.h>

namespace dials { namespace algorithms { namespace shoebox { namespace boost_python {

  using namespace boost::python;
  // testing
  void export_mask_builder() {
    def("build_mask",
        &build_mask,
        (arg("nx") = 23, arg("ny") = 17, arg("nrx") = 3, arg("nry") = 2, arg("nc") = 8),
        arg("data2d"));
  }

}}}}  // namespace dials::algorithms::shoebox::boost_python
