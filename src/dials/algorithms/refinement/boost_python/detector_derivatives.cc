// detector_derivatives.cc
// boost.python wrapper for compute_detector_derivatives.
// Registered into dials_refinement_helpers_ext.

#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/tuple.hpp>
#include <scitbx/array_family/boost_python/shared_flex_conversions.h>
#include "../parameterisation/detector_derivatives.h"

using namespace boost::python;

namespace af = scitbx::af;
using scitbx::mat3;
using scitbx::vec3;

namespace dials { namespace refinement { namespace boost_python {

  // Python-facing wrapper: converts af::shared<dXY_pair> to a boost::python::list
  // of 2-tuples (dX, dY) where each element is a flex.double array.
  //
  // dXY_pair is NOT exposed to Python as a class.  Using def_readonly on members
  // of type af::shared<double> fails with "No Python class registered for C++
  // class scitbx::af::shared<double>" because boost.python looks for a class_<>
  // wrapper rather than a to-python converter.  Returning tuples instead routes
  // through the to-python converter path, which works correctly once
  // shared_flex_conversions<double> has been instantiated below.
  boost::python::list compute_detector_derivatives_py(
    af::const_ref<mat3<double> > const& D,
    af::const_ref<vec3<double> > const& pv,
    af::const_ref<double> const& w_inv,
    af::const_ref<double> const& u_w_inv,
    af::const_ref<double> const& v_w_inv,
    boost::python::object const& dd_ddet_p) {
    af::shared<dXY_pair> result =
      compute_detector_derivatives(D, pv, w_inv, u_w_inv, v_w_inv, dd_ddet_p);
    boost::python::list py_result;
    for (std::size_t p = 0; p < result.size(); ++p) {
      // Each element is a (dX, dY) tuple of flex.double arrays.
      // pair[0] == dX, pair[1] == dY.  Empty arrays signal null parameters.
      py_result.append(boost::python::make_tuple(result[p].dX, result[p].dY));
    }
    return py_result;
  }

  void export_detector_derivatives() {
    // Register af::shared<double> <-> flex.double converters in this translation
    // unit.  The cctbx flex module registers them globally at import time, but
    // we instantiate here as well to be safe and self-contained.
    scitbx::af::boost_python::shared_flex_conversions<double>();

    def("compute_detector_derivatives",
        &compute_detector_derivatives_py,
        (arg("D"),
         arg("pv"),
         arg("w_inv"),
         arg("u_w_inv"),
         arg("v_w_inv"),
         arg("dd_ddet_p")));
  }

}}}  // namespace dials::refinement::boost_python
