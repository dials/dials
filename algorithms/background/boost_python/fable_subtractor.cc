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

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_fable_subtractor()
  {
    // Overloads for call method
    double (FableSubtractor::*call_shoebox_and_mask)(
      const af::const_ref< double, af::c_grid<3> >&, 
      af::ref< int, af::c_grid<3> >) const = &FableSubtractor::operator();
    void (FableSubtractor::*call_reflection)(Reflection &) const =
      &FableSubtractor::operator();  
    void (FableSubtractor::*call_reflection_list)(af::ref<Reflection>) const =
      &FableSubtractor::operator();  
  
    class_<FableSubtractor>("FableSubtractorAlgorithm", no_init)
      .def(init<std::size_t, double>((
        arg("min_data") = 10,
        arg("n_sigma") = 3.0)))
      .def("__call__", call_shoebox_and_mask, (
        arg("shoebox"),
        arg("mask")))
      .def("__call__", call_reflection, (
        arg("reflection")))
      .def("__call__", call_reflection_list, (
        arg("reflection_list")));
  }

}}} // namespace = dials::algorithms::boost_python
