#include <boost/python.hpp>
#include <boost/python/def.hpp>

using namespace boost::python;

namespace dials { namespace refinement { namespace boost_python {

  void export_prediction_parameter_helpers();

  BOOST_PYTHON_MODULE(dials_refinement_ext)
  {
    export_prediction_parameter_helpers();
  }

}}} // namespace dials::refinement::boost_python
