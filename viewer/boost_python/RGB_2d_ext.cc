#include <boost/python.hpp>
#include <boost/python/def.hpp>

namespace dials { namespace scratch { namespace boost_python {
  using namespace boost::python;

  void rgb_ext();

  BOOST_PYTHON_MODULE(rgb_ext)
  {
    rgb_ext();
  }

}}}
