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
#include <dials/algorithms/background/fable_subtractor.h>

namespace dials { namespace algorithms { namespace background {
  namespace boost_python {

  using namespace boost::python;

  static
  void subtract_with_shoebox(const FableSubtractor &subtract,
      Shoebox<> &shoebox) {

    typedef Shoebox<>::float_type FloatType;

    af::ref< FloatType, af::c_grid<3> > background = shoebox.background.ref();
    FloatType value = subtract(
      shoebox.data.const_ref(),
      shoebox.mask.ref());
    for (std::size_t i = 0; i < background.size(); ++i) {
      background[i] = value;
    }
  }

  static
  af::shared<bool> subtract_with_shoebox_list(
      const FableSubtractor &subtract,
      af::ref< Shoebox<> > shoeboxes) {
    af::shared<bool> result(shoeboxes.size(), true);
    for (std::size_t i = 0; i < shoeboxes.size(); ++i) {
      try {
        subtract_with_shoebox(subtract, shoeboxes[i]);
      } catch(dials::error) {
        result[i] = false;
      }
    }
    return result;
  }

  void export_fable_subtractor()
  {
    class_<FableSubtractor>("FableSubtractorAlgorithm", no_init)
      .def(init<std::size_t, double>((
        arg("min_data") = 10,
        arg("n_sigma") = 3.0)))
      .def("__call__",
        &FableSubtractor::operator()<float>, (
          arg("shoebox"),
          arg("mask")))
      .def("__call__",
        &FableSubtractor::operator()<double>, (
          arg("shoebox"),
          arg("mask")))
      .def("__call__", subtract_with_shoebox, (
        arg("shoebox")))
      .def("__call__", subtract_with_shoebox_list, (
        arg("shoebox_list")));
  }

}}}} // namespace = dials::algorithms::background::boost_python
