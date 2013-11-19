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
#include <dials/algorithms/integration/summation.h>
#include <dials/model/data/reflection.h>
#include <dials/algorithms/shoebox/mask_code.h>

using namespace boost::python;

namespace dials { namespace algorithms { namespace boost_python {

  using dials::model::Reflection;

  template <typename FloatType>
  void summation_wrapper(const char *name)
  {
    typedef Summation<FloatType> SummationType;

    class_ <SummationType> ("Summation", no_init)
      .def(init <const af::const_ref<FloatType>&,
                 const af::const_ref<FloatType>&>((
          arg("signal"),
          arg("background"))))
      .def(init <const af::const_ref<FloatType>&,
                 const af::const_ref<FloatType>&,
                 const af::const_ref<bool>&>((
          arg("signal"),
          arg("background"),
          arg("mask"))))
      .def(init <const af::const_ref< FloatType, af::c_grid<2> >&,
                 const af::const_ref< FloatType, af::c_grid<2> >&>((
          arg("signal"),
          arg("background"))))
      .def(init <const af::const_ref< FloatType, af::c_grid<2> >&,
                 const af::const_ref< FloatType, af::c_grid<2> >&,
                 const af::const_ref< bool, af::c_grid<2> >&>((
          arg("signal"),
          arg("background"),
          arg("mask"))))
      .def(init <const af::const_ref< FloatType, af::c_grid<3> >&,
                 const af::const_ref< FloatType, af::c_grid<3> >&>((
          arg("signal"),
          arg("background"))))
      .def(init <const af::const_ref< FloatType, af::c_grid<3> >&,
                 const af::const_ref< FloatType, af::c_grid<3> >&,
                 const af::const_ref< bool, af::c_grid<3> >&>((
          arg("signal"),
          arg("background"),
          arg("mask"))))
      .def("intensity",
        &SummationType::intensity)
      .def("variance",
        &SummationType::variance)
      .def("standard_deviation",
        &SummationType::standard_deviation)
      .def("signal_intensity",
        &SummationType::signal_intensity)
      .def("signal_variance",
        &SummationType::signal_variance)
      .def("signal_standard_deviation",
        &SummationType::signal_standard_deviation)
      .def("background_intensity",
        &SummationType::background_intensity)
      .def("background_variance",
        &SummationType::background_variance)
      .def("background_standard_deviation",
        &SummationType::background_standard_deviation);
  }

  template <typename FloatType>
  Summation<FloatType> make_summation_1d(
    const af::const_ref<FloatType> &image,
    const af::const_ref<FloatType> &background) {
    return Summation<FloatType>(image, background);
  }

  template <typename FloatType>
  Summation<FloatType> make_summation_1d_bg(
    const af::const_ref<FloatType> &image,
    const af::const_ref<FloatType> &background,
    const af::const_ref<bool> &mask) {
    return Summation<FloatType>(image, background, mask);
  }

  template <typename FloatType>
  Summation<FloatType> make_summation_2d(
    const af::const_ref<FloatType, af::c_grid<2> > &image,
    const af::const_ref<FloatType, af::c_grid<2> > &background) {
    return Summation<FloatType>(image, background);
  }

  template <typename FloatType>
  Summation<FloatType> make_summation_2d_bg(
    const af::const_ref<FloatType, af::c_grid<2> > &image,
    const af::const_ref<FloatType, af::c_grid<2> > &background,
    const af::const_ref<bool, af::c_grid<2> > &mask) {
    return Summation<FloatType>(image, background, mask);
  }

  template <typename FloatType>
  Summation<FloatType> make_summation_3d(
    const af::const_ref<FloatType, af::c_grid<3> > &image,
    const af::const_ref<FloatType, af::c_grid<3> > &background) {
    return Summation<FloatType>(image, background);
  }

  template <typename FloatType>
  Summation<FloatType> make_summation_3d_bg(
    const af::const_ref<FloatType, af::c_grid<3> > &image,
    const af::const_ref<FloatType, af::c_grid<3> > &background,
    const af::const_ref<bool, af::c_grid<3> > &mask) {
    return Summation<FloatType>(image, background, mask);
  }

  template <typename FloatType>
  void summation_suite() {

    def("integrate_by_summation",
      &make_summation_1d<FloatType>, (
        arg("image"),
        arg("background")));

    def("integrate_by_summation",
      &make_summation_1d_bg<FloatType>, (
        arg("image"),
        arg("background"),
        arg("mask")));

    def("integrate_by_summation",
      &make_summation_2d<FloatType>, (
        arg("image"),
        arg("background")));

    def("integrate_by_summation",
      &make_summation_2d_bg<FloatType>, (
        arg("image"),
        arg("background"),
        arg("mask")));

    def("integrate_by_summation",
      &make_summation_3d<FloatType>, (
        arg("image"),
        arg("background")));

    def("integrate_by_summation",
      &make_summation_3d_bg<FloatType>, (
        arg("image"),
        arg("background"),
        arg("mask")));
  }

  void summation_with_reflection(Reflection &r) {

    af::const_ref< int, af::c_grid<3> > shoebox_mask =
      r.get_shoebox_mask().const_ref();
    af::versa< bool, af::c_grid<3> > mask(shoebox_mask.accessor(),
      af::init_functor_null<bool>());
    for (std::size_t i = 0; i < mask.size(); ++i) {
      mask[i] = (shoebox_mask[i] & shoebox::Valid) ? true : false;
    }

    // Integrate the reflection
    Summation<Reflection::float_type> result(
      r.get_shoebox().const_ref(),
      r.get_shoebox_background().const_ref(),
      mask.const_ref());

    // Set the intensity and variance
    r.set_intensity(result.intensity());
    r.set_intensity_variance(result.variance());
  }

  void summation_with_reflection_list(af::ref<Reflection> reflections) {
    #pragma omp parallel for
    for (std::size_t i = 0; i < reflections.size(); ++i) {
      try {
        if (reflections[i].is_valid()) {
          summation_with_reflection(reflections[i]);
        }
      } catch (dials::error) {
        reflections[i].set_valid(false);
      }
    }
  }

  void export_summation() {
    summation_wrapper<float>("SummationFloat");
    summation_wrapper<double>("SummationDouble");

    summation_suite<float>();
    summation_suite<double>();

    def("integrate_by_summation", &summation_with_reflection);
    def("integrate_by_summation", &summation_with_reflection_list);
  }

}}} // namespace = dials::algorithms::boost_python
