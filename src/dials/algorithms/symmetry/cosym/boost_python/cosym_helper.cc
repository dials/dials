#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <cctbx/miller.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/algorithms/symmetry/cosym/cosym_helper.h>

namespace cosym { namespace boost_python {
  using namespace boost::python;

  void export_match_indices() {
    class_<matcher>("matcher", no_init)
      .def(init<scitbx::af::const_ref<cctbx::miller::index<> > const&>(
        (arg("miller_index"))))
      .def("match", &matcher::match);
  }
}}  // namespace cosym::boost_python