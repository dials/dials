#include <boost/python.hpp>
#include <boost/python/def.hpp>

namespace dials { namespace viewer { namespace boost_python {
  using namespace boost::python;
  void dials_viewer();

  BOOST_PYTHON_MODULE(dials_viewer_ext){
    dials_viewer();

  }

}}}
