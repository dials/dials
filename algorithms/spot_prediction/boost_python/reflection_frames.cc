/*
 * calculate_reflection_frames.cc
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
#include <dials/algorithms/spot_prediction/reflection_frames.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_reflection_frames()
  {
    // Function pointers to overloads
    shared<Reflection> (*reflection_frames_single)(const ScanData&,
      const Reflection&) = &reflection_frames;
    shared<Reflection> (*reflection_frames_array)(const ScanData&,
      const ReflectionList&) = &reflection_frames;

    // Export functions
    def("reflection_frames",
      reflection_frames_single, (
        arg("scan"),
        arg("reflection")));

    def("reflection_frames",
      reflection_frames_array, (
        arg("scan"),
        arg("reflection_list")));
  }

}}} // namespace = dials::algorithms::boost_python
