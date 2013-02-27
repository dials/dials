/*
 * spot_predictor.cc
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
#include <dials/algorithms/spot_prediction/spot_predictor.h>
#include <dxtbx/model/scan.h>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/goniometer.h>
#include <dxtbx/model/detector.h>
#include <dials/model/data/reflection.h>
#include "spot_predictor_wrapper.h"

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  using dxtbx::model::ScanData;
  using dxtbx::model::Beam;
  using dxtbx::model::Goniometer;
  using dxtbx::model::Detector;
  using model::Reflection;

  void export_spot_predictor()
  {
    // Export a spot predictor for flat panel detectors
    spot_predictor_wrapper <
      SpotPredictor <
        Beam,
        Detector,
        Goniometer,
        ScanData,
        Reflection> >("SpotPredictor");
  }

}}} // namespace = dials::spot_prediction::boost_python
