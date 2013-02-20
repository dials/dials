
#include <boost/python.hpp>
#include <boost/python/def.hpp>

namespace dials { namespace model { namespace boost_python {

  using namespace boost::python;

  void export_reflection();

  BOOST_PYTHON_MODULE(dials_reflection_ext)
  {
    export_reflection();
  }

}}} // namespace = dials::model::boost_python
