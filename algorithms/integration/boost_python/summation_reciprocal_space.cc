/*
 * summation_reciprocal_space.cc
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
#include <dials/algorithms/integration/summation_reciprocal_space.h>

using namespace boost::python;

namespace dials { namespace algorithms { namespace boost_python {

  void export_summation_reciprocal_space()
  {
    void (SummationReciprocalSpace::*call_w_reflection)(
      Reflection &r) const = &SummationReciprocalSpace::operator();
    void (SummationReciprocalSpace::*call_w_reflection_list)(
        af::ref<Reflection> reflections) const = 
          &SummationReciprocalSpace::operator();     
  
    class_<SummationReciprocalSpace>(
        "SummationReciprocalSpaceAlgorithm", no_init)
      .def(init<const shared_ptr<Beam>&,
                const shared_ptr<Detector>&,
                const shared_ptr<Goniometer>&>((
        arg("beam"),
        arg("detector"),
        arg("goniometer"))))
      .def("__call__", call_w_reflection)
      .def("__call__", call_w_reflection_list);
  }

}}} // namespace = dials::algorithms::boost_python
