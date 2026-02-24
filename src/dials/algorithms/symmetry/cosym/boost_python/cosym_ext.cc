#include <boost/python.hpp>
#include <boost/python/def.hpp>

using namespace boost::python;

namespace cosym { namespace boost_python {
  void export_match_indices();

  BOOST_PYTHON_MODULE(dials_cosym_ext) {
    export_match_indices();
  }

}}  // namespace cosym::boost_python