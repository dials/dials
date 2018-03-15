#include <boost/python.hpp>
#include <boost/python/def.hpp>

using namespace boost::python;

namespace dials_scratch { namespace scaling { namespace boost_python {

  void export_elementwise_square();
  void export_sph_harm_table();

  BOOST_PYTHON_MODULE(dials_scratch_scaling_ext)
  {
    export_elementwise_square();
    export_sph_harm_table();
  }

}}} // namespace dials_scratch::scaling::boost_python
