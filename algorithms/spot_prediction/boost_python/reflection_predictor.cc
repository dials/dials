/*
 * reflection_predictor.cc
 *
 *  Copyright (C) 2014 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/algorithms/spot_prediction/reflection_predictor.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_scan_static_reflection_predictor() {

    af::reflection_table (ScanStaticReflectionPredictor::*predict_all)()const = &ScanStaticReflectionPredictor::operator();
    af::reflection_table (ScanStaticReflectionPredictor::*predict_observed)(const af::const_ref< cctbx::miller::index<> >&)const = &ScanStaticReflectionPredictor::operator();

    af::reflection_table (ScanStaticReflectionPredictor::*predict_observed_with_panel)(const af::const_ref< cctbx::miller::index<> >&, std::size_t)const = &ScanStaticReflectionPredictor::operator();

    af::reflection_table (ScanStaticReflectionPredictor::*predict_observed_with_panel_list)(const af::const_ref< cctbx::miller::index<> >&, const af::const_ref<std::size_t>&)const = &ScanStaticReflectionPredictor::operator();

    class_<ScanStaticReflectionPredictor>("ScanStaticReflectionPredictor", no_init)
      .def(init<
          const Beam&,
          const Detector&,
          const Goniometer&,
          const Scan&,
          const cctbx::uctbx::unit_cell&,
          const cctbx::sgtbx::space_group_type&,
          mat3<double>,
          double>())
      .def("__call__", predict_all)
      .def("__call__", predict_observed)
      .def("__call__", predict_observed_with_panel)
      .def("__call__", predict_observed_with_panel_list);
  }

  void export_scan_varying_reflection_predictor() {


    af::reflection_table (ScanVaryingReflectionPredictor::*predict_all)()const = &ScanVaryingReflectionPredictor::operator();
    af::reflection_table (ScanVaryingReflectionPredictor::*predict_observed)(const af::const_ref< cctbx::miller::index<> >&, const af::const_ref<int>&)const = &ScanVaryingReflectionPredictor::operator();

    af::reflection_table (ScanVaryingReflectionPredictor::*predict_observed_with_panel)(const af::const_ref< cctbx::miller::index<> >&, const af::const_ref<int>&, std::size_t)const = &ScanVaryingReflectionPredictor::operator();

    af::reflection_table (ScanVaryingReflectionPredictor::*predict_observed_with_panel_list)(const af::const_ref< cctbx::miller::index<> >&, const af::const_ref<int>&, const af::const_ref<std::size_t>&)const = &ScanVaryingReflectionPredictor::operator();

    class_<ScanVaryingReflectionPredictor>("ScanVaryingReflectionPredictor", no_init)
      .def(init<
          const Beam&,
          const Detector&,
          const Goniometer&,
          const Scan&,
          const af::const_ref< mat3<double> > &,
          double,
          std::size_t>())
      .def("__call__", predict_all)
      .def("__call__", predict_observed)
      .def("__call__", predict_observed_with_panel)
      .def("__call__", predict_observed_with_panel_list)
      ;
  }

  void export_reflection_predictor()
  {
    export_scan_static_reflection_predictor();
    export_scan_varying_reflection_predictor();

    class_<StillsReflectionPredictor>("StillsReflectionPredictor", no_init)
      .def("observed", &StillsReflectionPredictor::observed);
  }

}}} // namespace = dials::algorithms::boost_python
