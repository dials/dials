/*
 * xds_subtractor.cc
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
#include <dials/algorithms/background/xds_subtractor.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  static
  void subtract_with_reflection(const XdsSubtractor &subtract,
      Reflection &reflection) {

    typedef Reflection::float_type FloatType;

    af::ref< int, af::c_grid<3> > mask = reflection.get_shoebox_mask().ref();
    af::ref< FloatType, af::c_grid<3> > shoebox = reflection.get_shoebox().ref();
    af::ref< FloatType, af::c_grid<3> > background =
      reflection.get_shoebox_background().ref();

    FloatType value = subtract(shoebox, mask);
    for (std::size_t i = 0; i < background.size(); ++i) {
      background[i] = value;
    }
  }

  static
  void subtract_with_reflection_list(const XdsSubtractor &subtract,
      af::ref<Reflection> reflections) {
    #pragma omp parallel for
    for (std::size_t i = 0; i < reflections.size(); ++i) {
      try {
        if (reflections[i].is_valid()) {
          subtract_with_reflection(subtract, reflections[i]);
        }
      } catch(dials::error) {
        reflections[i].set_valid(false);
      }
    }
  }

  void export_xds_subtractor()
  {
    class_<XdsSubtractor>("XdsSubtractorAlgorithm", no_init)
      .def(init<std::size_t, double>((
        arg("min_data") = 10,
        arg("n_sigma") = 3.0)))
      .def("__call__",
         &XdsSubtractor::operator()<float>, (
          arg("shoebox"),
          arg("mask")))
      .def("__call__",
         &XdsSubtractor::operator()<double>, (
          arg("shoebox"),
          arg("mask")))
      .def("__call__",
        subtract_with_reflection, (
          arg("reflection")))
      .def("__call__",
        subtract_with_reflection_list, (
          arg("reflection_list")));
  }

}}} // namespace = dials::algorithms::boost_python
