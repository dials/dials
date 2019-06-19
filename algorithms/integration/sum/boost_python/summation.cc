/*
 * summation.cc
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
#include <dials/algorithms/integration/sum/summation.h>
#include <dials/algorithms/integration/sum/sum_image_volume.h>
#include <dials/algorithms/shoebox/mask_code.h>

using namespace boost::python;

namespace dials { namespace algorithms { namespace boost_python {

  template <typename FloatType>
  void summation_wrapper(const char *name) {
    typedef Summation<FloatType> SummationType;

    class_<SummationType>(name, no_init)
      .def(init<const af::const_ref<FloatType> &,
                const af::const_ref<FloatType> &,
                const af::const_ref<int> &>((boost::python::arg("signal"),
                                             boost::python::arg("background"),
                                             boost::python::arg("mask"))))
      .def(init<const af::const_ref<FloatType, af::c_grid<2> > &,
                const af::const_ref<FloatType, af::c_grid<2> > &,
                const af::const_ref<int, af::c_grid<2> > &>(
        (boost::python::arg("signal"),
         boost::python::arg("background"),
         boost::python::arg("mask"))))
      .def(init<const af::const_ref<FloatType, af::c_grid<3> > &,
                const af::const_ref<FloatType, af::c_grid<3> > &,
                const af::const_ref<int, af::c_grid<3> > &>(
        (boost::python::arg("signal"),
         boost::python::arg("background"),
         boost::python::arg("mask"))))
      .def("intensity", &SummationType::intensity)
      .def("variance", &SummationType::variance)
      .def("background", &SummationType::background)
      .def("background_variance", &SummationType::background_variance)
      .def("n_signal", &SummationType::n_signal)
      .def("n_background", &SummationType::n_background)
      .def("success", &SummationType::success);
  }

  template <typename FloatType>
  Summation<FloatType> make_summation_1d(const af::const_ref<FloatType> &image,
                                         const af::const_ref<FloatType> &background,
                                         const af::const_ref<int> &mask) {
    return Summation<FloatType>(image, background, mask);
  }

  template <typename FloatType>
  Summation<FloatType> make_summation_2d(
    const af::const_ref<FloatType, af::c_grid<2> > &image,
    const af::const_ref<FloatType, af::c_grid<2> > &background,
    const af::const_ref<int, af::c_grid<2> > &mask) {
    return Summation<FloatType>(image, background, mask);
  }

  template <typename FloatType>
  Summation<FloatType> make_summation_3d(
    const af::const_ref<FloatType, af::c_grid<3> > &image,
    const af::const_ref<FloatType, af::c_grid<3> > &background,
    const af::const_ref<int, af::c_grid<3> > &mask) {
    return Summation<FloatType>(image, background, mask);
  }

  template <typename FloatType>
  void summation_suite() {
    def("integrate_by_summation",
        &make_summation_1d<FloatType>,
        (boost::python::arg("image"),
         boost::python::arg("background"),
         boost::python::arg("mask")));

    def("integrate_by_summation",
        &make_summation_2d<FloatType>,
        (boost::python::arg("image"),
         boost::python::arg("background"),
         boost::python::arg("mask")));

    def("integrate_by_summation",
        &make_summation_3d<FloatType>,
        (boost::python::arg("image"),
         boost::python::arg("background"),
         boost::python::arg("mask")));
  }

  void export_summation() {
    summation_wrapper<float>("SummationFloat");
    summation_wrapper<double>("SummationDouble");

    summation_suite<float>();
    summation_suite<double>();

    def("sum_image_volume", &sum_multi_panel_image_volume<float>);
  }

}}}  // namespace dials::algorithms::boost_python
