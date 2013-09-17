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

using namespace boost::python;

namespace dials { namespace algorithms { namespace boost_python {

  void export_summation()
  {
    class_ <Summation> ("Summation", no_init)
      .def(init <const af::const_ref<double>&,
                 const af::const_ref<double>&>((
          arg("signal"),
          arg("background"))))
      .def(init <const af::const_ref<double>&,
                 const af::const_ref<double>&,
                 const af::const_ref<bool>&>((
          arg("signal"),
          arg("background"),
          arg("mask"))))
      .def(init <const af::const_ref< double, af::c_grid<2> >&,
                 const af::const_ref< double, af::c_grid<2> >&>((
          arg("signal"),
          arg("background"))))
      .def(init <const af::const_ref< double, af::c_grid<2> >&,
                 const af::const_ref< double, af::c_grid<2> >&,
                 const af::const_ref< bool, af::c_grid<2> >&>((
          arg("signal"),
          arg("background"),
          arg("mask"))))
      .def(init <const af::const_ref< double, af::c_grid<3> >&,
                 const af::const_ref< double, af::c_grid<3> >&>((
          arg("signal"),
          arg("background"))))
      .def(init <const af::const_ref< double, af::c_grid<3> >&,
                 const af::const_ref< double, af::c_grid<3> >&,
                 const af::const_ref< bool, af::c_grid<3> >&>((
          arg("signal"),
          arg("background"),
          arg("mask"))))
      .def("intensity", &Summation::intensity)
      .def("variance", &Summation::variance)
      .def("standard_deviation", &Summation::standard_deviation)
      .def("signal_intensity", &Summation::signal_intensity)
      .def("signal_variance", &Summation::signal_variance)
      .def("signal_standard_deviation", &Summation::signal_standard_deviation)
      .def("background_intensity", &Summation::background_intensity)
      .def("background_variance", &Summation::background_variance)
      .def("background_standard_deviation", 
        &Summation::background_standard_deviation);
        
    Summation3d::integrator (Summation3d::*call_w_pixels)(
        const af::const_ref<double, af::c_grid<3> > &pixels, 
        const af::const_ref<double, af::c_grid<3> > &background,
        const af::const_ref<bool, af::c_grid<3> > &mask) const = 
      &Summation3d::operator();
    void (Summation3d::*call_w_reflection)(Reflection &r) const = 
      &Summation3d::operator();
    void (Summation3d::*call_w_reflection_list)(
        af::ref<Reflection> reflections) const = 
      &Summation3d::operator();         
        
    class_<Summation3d>("Summation3dAlgorithm")
      .def("__call__", call_w_pixels)
      .def("__call__", call_w_reflection)
      .def("__call__", call_w_reflection_list);
  }

}}} // namespace = dials::algorithms::boost_python
