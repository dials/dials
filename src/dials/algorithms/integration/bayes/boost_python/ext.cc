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
#include <dials/algorithms/integration/bayes/bayesian_integrator.h>
#include <dials/algorithms/shoebox/mask_code.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  template <typename FloatType>
  void bayesian_integrator_wrapper(const char *name) {
    typedef BayesianIntegrator<FloatType> BayesianIntegratorType;

    class_<BayesianIntegratorType>(name, no_init)
      .def(init<const af::const_ref<FloatType> &,
                const af::const_ref<FloatType> &,
                const af::const_ref<int> &>(
        (arg("signal"), arg("background"), arg("mask"))))
      .def(init<const af::const_ref<FloatType, af::c_grid<2> > &,
                const af::const_ref<FloatType, af::c_grid<2> > &,
                const af::const_ref<int, af::c_grid<2> > &>(
        (arg("signal"), arg("background"), arg("mask"))))
      .def(init<const af::const_ref<FloatType, af::c_grid<3> > &,
                const af::const_ref<FloatType, af::c_grid<3> > &,
                const af::const_ref<int, af::c_grid<3> > &>(
        (arg("signal"), arg("background"), arg("mask"))))
      .def("intensity", &BayesianIntegratorType::intensity)
      .def("variance", &BayesianIntegratorType::variance)
      .def("background", &BayesianIntegratorType::background)
      .def("background_variance", &BayesianIntegratorType::background_variance)
      .def("n_signal", &BayesianIntegratorType::n_signal)
      .def("n_background", &BayesianIntegratorType::n_background)
      .def("success", &BayesianIntegratorType::success);
  }

  template <typename FloatType>
  BayesianIntegrator<FloatType> make_bayesian_integrator_1d(
    const af::const_ref<FloatType> &image,
    const af::const_ref<FloatType> &background,
    const af::const_ref<int> &mask) {
    return BayesianIntegrator<FloatType>(image, background, mask);
  }

  template <typename FloatType>
  BayesianIntegrator<FloatType> make_bayesian_integrator_2d(
    const af::const_ref<FloatType, af::c_grid<2> > &image,
    const af::const_ref<FloatType, af::c_grid<2> > &background,
    const af::const_ref<int, af::c_grid<2> > &mask) {
    return BayesianIntegrator<FloatType>(image, background, mask);
  }

  template <typename FloatType>
  BayesianIntegrator<FloatType> make_bayesian_integrator_3d(
    const af::const_ref<FloatType, af::c_grid<3> > &image,
    const af::const_ref<FloatType, af::c_grid<3> > &background,
    const af::const_ref<int, af::c_grid<3> > &mask) {
    return BayesianIntegrator<FloatType>(image, background, mask);
  }

  template <typename FloatType>
  void bayesian_integrator_suite() {
    def("integrate_by_bayesian_integrator",
        &make_bayesian_integrator_1d<FloatType>,
        (arg("image"), arg("background"), arg("mask")));

    def("integrate_by_bayesian_integrator",
        &make_bayesian_integrator_2d<FloatType>,
        (arg("image"), arg("background"), arg("mask")));

    def("integrate_by_bayesian_integrator",
        &make_bayesian_integrator_3d<FloatType>,
        (arg("image"), arg("background"), arg("mask")));
  }

  BOOST_PYTHON_MODULE(dials_algorithms_integration_bayes_ext) {
    bayesian_integrator_wrapper<float>("BayesianIntegratorFloat");
    bayesian_integrator_wrapper<double>("BayesianIntegratorDouble");

    bayesian_integrator_suite<float>();
    bayesian_integrator_suite<double>();
  }

}}}  // namespace dials::algorithms::boost_python
