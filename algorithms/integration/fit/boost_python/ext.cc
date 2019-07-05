/*
 * fitting.cc
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
#include <dials/algorithms/integration/fit/fitting.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  template <typename FloatType>
  void profile_fitter_wrapper(const char *name) {
    typedef ProfileFitter<FloatType> ProfileFitterType;

    class_<ProfileFitterType>(name, no_init)
      .def("intensity", &ProfileFitterType::intensity)
      .def("variance", &ProfileFitterType::variance)
      .def("correlation", &ProfileFitterType::correlation)
      .def("niter", &ProfileFitterType::niter)
      .def("maxiter", &ProfileFitterType::maxiter)
      .def("error", &ProfileFitterType::error);
    ;
  }

  template <typename FloatType>
  ProfileFitter<FloatType> make_profile_fitter_1d_1(const af::const_ref<FloatType> &d,
                                                    const af::const_ref<FloatType> &b,
                                                    const af::const_ref<bool> &m,
                                                    const af::const_ref<FloatType> &p,
                                                    double eps,
                                                    std::size_t maxiter) {
    return ProfileFitter<FloatType>(d, b, m, p, eps, maxiter);
  }

  template <typename FloatType>
  ProfileFitter<FloatType> make_profile_fitter_2d_1(
    const af::const_ref<FloatType, af::c_grid<2> > &d,
    const af::const_ref<FloatType, af::c_grid<2> > &b,
    const af::const_ref<bool, af::c_grid<2> > &m,
    const af::const_ref<FloatType, af::c_grid<2> > &p,
    double eps,
    std::size_t maxiter) {
    return ProfileFitter<FloatType>(d, b, m, p, eps, maxiter);
  }

  template <typename FloatType>
  ProfileFitter<FloatType> make_profile_fitter_3d_1(
    const af::const_ref<FloatType, af::c_grid<3> > &d,
    const af::const_ref<FloatType, af::c_grid<3> > &b,
    const af::const_ref<bool, af::c_grid<3> > &m,
    const af::const_ref<FloatType, af::c_grid<3> > &p,
    double eps,
    std::size_t maxiter) {
    return ProfileFitter<FloatType>(d, b, m, p, eps, maxiter);
  }

  template <typename FloatType>
  ProfileFitter<FloatType> make_profile_fitter_1d_n(
    const af::const_ref<FloatType> &d,
    const af::const_ref<FloatType> &b,
    const af::const_ref<bool> &m,
    const af::const_ref<FloatType, af::c_grid<2> > &p,
    double eps,
    std::size_t maxiter) {
    return ProfileFitter<FloatType>(d, b, m, p, eps, maxiter);
  }

  template <typename FloatType>
  ProfileFitter<FloatType> make_profile_fitter_2d_n(
    const af::const_ref<FloatType, af::c_grid<2> > &d,
    const af::const_ref<FloatType, af::c_grid<2> > &b,
    const af::const_ref<bool, af::c_grid<2> > &m,
    const af::const_ref<FloatType, af::c_grid<3> > &p,
    double eps,
    std::size_t maxiter) {
    return ProfileFitter<FloatType>(d, b, m, p, eps, maxiter);
  }

  template <typename FloatType>
  ProfileFitter<FloatType> make_profile_fitter_3d_n(
    const af::const_ref<FloatType, af::c_grid<3> > &d,
    const af::const_ref<FloatType, af::c_grid<3> > &b,
    const af::const_ref<bool, af::c_grid<3> > &m,
    const af::const_ref<FloatType, af::c_grid<4> > &p,
    double eps,
    std::size_t maxiter) {
    return ProfileFitter<FloatType>(d, b, m, p, eps, maxiter);
  }

  template <typename Func>
  void def_make_profile_fitter(Func func) {
    def("ProfileFitter",
        func,
        (arg("data"),
         arg("background"),
         arg("mask"),
         arg("profile"),
         arg("eps") = 1e-3,
         arg("maxiter") = 10));
  }

  BOOST_PYTHON_MODULE(dials_algorithms_integration_fit_ext) {
    profile_fitter_wrapper<float>("ProfileFitterFloat");
    profile_fitter_wrapper<double>("ProfileFitterDouble");

    def_make_profile_fitter(&make_profile_fitter_1d_1<float>);
    def_make_profile_fitter(&make_profile_fitter_2d_1<float>);
    def_make_profile_fitter(&make_profile_fitter_2d_1<float>);
    def_make_profile_fitter(&make_profile_fitter_1d_n<float>);
    def_make_profile_fitter(&make_profile_fitter_2d_n<float>);
    def_make_profile_fitter(&make_profile_fitter_3d_n<float>);

    def_make_profile_fitter(&make_profile_fitter_1d_1<double>);
    def_make_profile_fitter(&make_profile_fitter_2d_1<double>);
    def_make_profile_fitter(&make_profile_fitter_3d_1<double>);
    def_make_profile_fitter(&make_profile_fitter_1d_n<double>);
    def_make_profile_fitter(&make_profile_fitter_2d_n<double>);
    def_make_profile_fitter(&make_profile_fitter_3d_n<double>);
  }

}}}  // namespace dials::algorithms::boost_python
