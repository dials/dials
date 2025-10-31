
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/algorithms/integration/tof/tof_mask_calculator.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;
  BOOST_PYTHON_MODULE(dials_algorithms_tof_integration_ext) {
    def("tof_calculate_ellipse_shoebox_mask",
        &tof_calculate_ellipse_shoebox_mask,
        (arg("reflection_table"), arg("experiment")));

    def("tof_calculate_seed_skewness_shoebox_mask",
        &tof_calculate_seed_skewness_shoebox_mask,
        (arg("reflection_table"),
         arg("experiment"),
         arg("d_skewness_threshold"),
         arg("min_iterations")));
  }

}}}  // namespace dials::algorithms::boost_python
