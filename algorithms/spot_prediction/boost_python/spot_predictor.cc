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

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_spot_predictor()
  {
    scitbx::af::shared <Reflection> (SpotPredictor::*predict_single)(
      miller_index) const = &SpotPredictor::operator();
    scitbx::af::shared <Reflection> (SpotPredictor::*predict_array)(
      const flex_miller_index &) const = &SpotPredictor::operator();
    scitbx::af::shared <Reflection> (SpotPredictor::*predict_generate)() = 
      &SpotPredictor::operator();
                
    class_ <SpotPredictor> ("SpotPredictor", no_init)
      .def(init <const Beam &,
                 const FlatPanelDetector &,
                 const Goniometer &,
                 const Scan &,
                 const cctbx::uctbx::unit_cell &,
                 const cctbx::sgtbx::space_group_type &,
                 mat3 <double>,
                 double> ((
        arg("beam"),
        arg("detector"),
        arg("goniometer"),
        arg("scan"),
        arg("unit_cell"),
        arg("space_group_type"),
        arg("ub_matrix"),
        arg("d_min"))))
      .def("__call__", predict_single, (
        arg("miller_index")))
      .def("__call__", predict_array, (
        arg("miller_indices")))
      .def("__call__", predict_generate);
  }

}}} // namespace = dials::spot_prediction::boost_python
