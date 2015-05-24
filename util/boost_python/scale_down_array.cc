/*
 * FIXME add a header
 */

#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/util/scale_down_array.h>

namespace dials { namespace util { namespace boost_python {
  using namespace boost::python;
  BOOST_PYTHON_MODULE(dials_util_ext) {
    def("scale_down_array", &scale_down_array,
        (arg("image"), arg("scale_factor")));
  }
}}}
