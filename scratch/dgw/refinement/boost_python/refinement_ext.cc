#include <boost/python.hpp>
#include <boost/python/def.hpp>

using namespace boost::python;

namespace bpcx_regression { namespace refinement { namespace boost_python {

  void export_prediction_parameter_helpers();

  BOOST_PYTHON_MODULE(bpcx_regression_refinement_ext)
  {
    export_prediction_parameter_helpers();
  }

}}} // namespace bpcx_regression::refinement::boost_python
