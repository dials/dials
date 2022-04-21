/*
 * stills_ray_predictor.cc
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
#include <dials/algorithms/spot_prediction/stills_ray_predictor.h>
#include <dxtbx/model/scan.h>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/goniometer.h>
#include <dxtbx/model/detector.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_stills_ray_predictor() {
    // Create and return the wrapper for the spot predictor object
    class_<StillsRayPredictor>("StillsRayPredictor", no_init)
      .def(init<vec3<double> >((arg("s0"))))
      .def(
        "__call__", &StillsRayPredictor::operator(), (arg("miller_index"), arg("UB")));
  }

}}}  // namespace dials::algorithms::boost_python
