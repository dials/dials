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
    // Useful typedefs
    typedef SpotPredictor::beam_type beam_type;
    typedef SpotPredictor::detector_type detector_type;
    typedef SpotPredictor::goniometer_type goniometer_type;
    typedef SpotPredictor::scan_type scan_type;
    typedef SpotPredictor::reflection_type reflection_type;
    typedef SpotPredictor::reflection_list_type reflection_list_type;

    // Typedef the different overloads for operator()
    reflection_list_type (SpotPredictor::*predict_single)(
      miller_index) const = &SpotPredictor::operator();
    reflection_list_type (SpotPredictor::*predict_array)(
      const flex_miller_index &) const = &SpotPredictor::operator();
    reflection_list_type (SpotPredictor::*predict_generate)() =
      &SpotPredictor::operator();

    // Create and return the wrapper for the spot predictor object
    class_ <SpotPredictor> ("SpotPredictor", no_init)
      .def(init <const beam_type&,
                 const detector_type&,
                 const goniometer_type&,
                 const scan_type&,
                 const cctbx::uctbx::unit_cell &,
                 const cctbx::sgtbx::space_group_type &,
                 mat3 <double>,
                 double> ((
        arg("beam"),
        arg("detector"),
        arg("goniometer"),
        arg("scan"),
        arg("unit_cell"),
        arg("space_group"),
        arg("UB"),
        arg("d_min"))))
      .def("__call__", predict_single, (
        arg("miller_index")))
      .def("__call__", predict_array, (
        arg("miller_indices")))
      .def("__call__", predict_generate);
  }

}}} // namespace = dials::spot_prediction::boost_python
