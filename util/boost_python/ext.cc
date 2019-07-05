/*
 * FIXME add a header
 */

#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/util/scale_down_array.h>
#include <dials/util/masking.h>
#include <dials/util/export_mtz_helpers.h>

namespace dials { namespace util { namespace boost_python {

  using namespace boost::python;
  BOOST_PYTHON_MODULE(dials_util_ext) {
    def("scale_down_array", &scale_down_array, (arg("image"), arg("scale_factor")));

    def("dials_u_to_mosflm", &dials_u_to_mosflm, (arg("dials_U"), arg("uc")));

    class_<ResolutionMaskGenerator>("ResolutionMaskGenerator", no_init)
      .def(init<const BeamBase &, const Panel &>())
      .def("apply", &ResolutionMaskGenerator::apply);
  }
}}}  // namespace dials::util::boost_python
