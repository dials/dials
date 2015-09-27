/*
 * reflection_predictor.cc
 *
 *  Copyright (C) 2015 Diamond Light Source
 *
 *  Author: Nicholas Sauter
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/scratch/idy/algorithms/spot_prediction/reflection_predictor.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_experimental_reflection_predictor() {

    typedef StillsExperimentalReflectionPredictor Predictor;

    af::reflection_table (Predictor::*predict_all)() const =
      &Predictor::operator();

    af::reflection_table (Predictor::*predict_observed)(
        const af::const_ref< cctbx::miller::index<> >&) =
      &Predictor::operator();

    af::reflection_table (Predictor::*predict_observed_with_panel)(
        const af::const_ref< cctbx::miller::index<> >&,
        std::size_t) = &Predictor::operator();

    af::reflection_table (Predictor::*predict_observed_with_panel_list)(
        const af::const_ref< cctbx::miller::index<> >&,
        const af::const_ref<std::size_t>&) = &Predictor::operator();

    class_<Predictor>("StillsExperimentalReflectionPredictor", no_init)
      .def(init<
          const Beam&,
          const Detector&,
          mat3<double>,
          const cctbx::uctbx::unit_cell&,
          const cctbx::sgtbx::space_group_type&,
          const double&,double const&>())
      .def("__call__", predict_all)
      .def("for_ub", &Predictor::for_ub)
      .def("__call__", predict_observed)
      .def("__call__", predict_observed_with_panel)
      .def("__call__", predict_observed_with_panel_list)
      .def("for_reflection_table",
          &Predictor::for_reflection_table)
      ;
  }

}}} // namespace = dials::algorithms::boost_python
