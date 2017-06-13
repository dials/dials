/*
 * ext.cc
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/iterator.hpp>
#include <dials/algorithms/integration/maximum_likelihood/fitting.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  template <typename FloatType>
  void profile_fitting_wrapper(const char *name) {

    typedef MLProfileFitting<FloatType> ProfileFittingType;

    class_<ProfileFittingType>(name, no_init)
      .def(init<const af::const_ref<FloatType, af::c_grid<3> >&,
                const af::const_ref<bool, af::c_grid<3> >&,
                const af::const_ref<FloatType, af::c_grid<3> >&,
                const af::const_ref<FloatType, af::c_grid<3> >&,
                double,
                std::size_t>((
        arg("profile"),
        arg("mask"),
        arg("contents"),
        arg("background"),
        arg("bits") = 1e-3,
        arg("max_iter") = 10)))
      .def("intensity", &ProfileFittingType::intensity)
      .def("variance", &ProfileFittingType::variance)
      .def("correlation", &ProfileFittingType::correlation)
      .def("niter", &ProfileFittingType::niter)
      .def("error", &ProfileFittingType::error);
  }

  template <typename FloatType>
  MLProfileFitting<FloatType> make_profile_fitting(
      const af::const_ref<FloatType, af::c_grid<3> > &p,
      const af::const_ref<bool, af::c_grid<3> > &m,
      const af::const_ref<FloatType, af::c_grid<3> > &c,
      const af::const_ref<FloatType, af::c_grid<3> > &b,
      double eps,
      std::size_t max_iter) {
    return MLProfileFitting<FloatType>(p, m, c, b, eps, max_iter);
  }

  template <typename FloatType>
  void profile_fitting_suite() {
    def("fit_profile", &make_profile_fitting<FloatType>, (
      arg("profile"),
      arg("mask"),
      arg("contents"),
      arg("background"),
      arg("bits") = 1e-3,
      arg("max_iter") = 10));
  }

  BOOST_PYTHON_MODULE(dials_algorithms_integration_maximum_likelihood_ext)
  {
    profile_fitting_wrapper<float>("ProfileFittingFloat");
    profile_fitting_wrapper<double>("ProfileFittingDouble");

    profile_fitting_suite<float>();
    profile_fitting_suite<double>();
  }

}}} // namespace = dials::algorithms::boost_python
