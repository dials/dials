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
#include <dials/algorithms/spot_prediction/multi_panel_spot_predictor.h>
#include <dials/model/data/reflection.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_multi_panel_spot_predictor()
  {
    // Useful typedefs
    typedef MultiPanelSpotPredictor::beam_type beam_type;
    typedef MultiPanelSpotPredictor::detector_type detector_type;
    typedef MultiPanelSpotPredictor::goniometer_type goniometer_type;
    typedef MultiPanelSpotPredictor::scan_type scan_type;
    typedef MultiPanelSpotPredictor::reflection_type reflection_type;
    typedef MultiPanelSpotPredictor::reflection_list_type reflection_list_type;

    // Typedef the different overloads for operator()
    reflection_list_type (MultiPanelSpotPredictor::*predict_single)(
      miller_index) const = &MultiPanelSpotPredictor::operator();
    reflection_list_type (MultiPanelSpotPredictor::*predict_array)(
      const flex_miller_index &) const = &MultiPanelSpotPredictor::operator();

    // Create and return the wrapper for the spot predictor object
    class_ <MultiPanelSpotPredictor> ("MultiPanelSpotPredictor", no_init)
      .def(init <const beam_type&,
                 const detector_type&,
                 const goniometer_type&,
                 const scan_type&,
                 mat3 <double> > ((
        arg("beam"),
        arg("detector"),
        arg("goniometer"),
        arg("scan"),
        arg("UB"))))
      .def("__call__", predict_single, (
        arg("miller_index")))
      .def("__call__", predict_array, (
        arg("miller_indices")));
  }

}}} // namespace = dials::spot_prediction::boost_python
