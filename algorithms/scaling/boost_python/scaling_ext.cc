#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/error.h>

using namespace boost::python;

namespace dials { namespace scaling { namespace boost_python {

  void export_grouped_obs();
  void export_minimum_multiplicity_selection();

  BOOST_PYTHON_MODULE(dials_scaling_helpers_ext)
  {
    export_grouped_obs();
    export_minimum_multiplicity_selection();
  }

}}} // namespace dials::refinement::boost_python
