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
    class_ <SumIntensity3d> ("SumIntensity3d", no_init)
      .def(init <const flex_double&,
                 const flex_double&>((
          arg("signal"),
          arg("background"))))
      .def(init <const flex_double&,
                 const flex_double&,
                 const flex_bool&>((
          arg("signal"),
          arg("background"),
          arg("mask"))))
      .def("intensity", &SumIntensity3d::intensity)
      .def("variance", &SumIntensity3d::variance)
      .def("signal_intensity", &SumIntensity3d::signal_intensity)
      .def("signal_variance", &SumIntensity3d::signal_variance)
      .def("background_intensity", &SumIntensity3d::background_intensity)
      .def("background_variance", &SumIntensity3d::background_variance);
      
    class_ <IntegrateBySummation> ("IntegrateBySummation", no_init)
      .def(init <const flex_double&, 
                 const flex_double&> ((
          arg("pixels"),
          arg("background"))))
      .def(init <const flex_double&, 
                 const flex_double&,
                 const flex_bool&> ((
          arg("pixels"),
          arg("mask"),
          arg("background"))))
      .def("intensity", 
        &IntegrateBySummation::intensity)
      .def("variance",
        &IntegrateBySummation::variance)
      .def("standard_deviation",
        &IntegrateBySummation::standard_deviation);
        
    Summation3d::integrator (Summation3d::*call_w_pixels)(
        const flex_double &pixels, 
        const flex_double &background,
        const flex_bool &mask) const = 
      &Summation3d::operator();
    void (Summation3d::*call_w_reflection)(Reflection &r) const = 
      &Summation3d::operator();
    void (Summation3d::*call_w_reflection_list)(
        ReflectionList &reflections) const = 
      &Summation3d::operator();         
        
    class_<Summation3d>("Summation3dAlgorithm")
      .def("__call__", call_w_pixels)
      .def("__call__", call_w_reflection)
      .def("__call__", call_w_reflection_list);
  }

}}} // namespace = dials::algorithms::boost_python
