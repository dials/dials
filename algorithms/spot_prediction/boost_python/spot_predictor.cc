/*
 * scan_predictor.cc
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
#include "spot_predictor_wrapper.h"

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_spot_predictor()
  {
    spot_predictor_wrapper <
      SpotPredictor <
        FlatPanelDetector, 
        Reflection> >("SpotPredictor");
            
    spot_predictor_wrapper <
      SpotPredictor <
        MultiFlatPanelDetector, 
        MultiPanelDetectorReflection> >("SpotPredictor2");  
  }

}}} // namespace = dials::spot_prediction::boost_python
