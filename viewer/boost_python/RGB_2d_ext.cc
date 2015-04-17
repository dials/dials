#include <boost/python.hpp>
#include <boost/python/def.hpp>

namespace dials { namespace scratch { namespace boost_python {
  using namespace boost::python;

  void dials_viewer_ext();

  BOOST_PYTHON_MODULE(dials_viewer_ext)
  {
    dials_viewer_ext();
  }

}}}
