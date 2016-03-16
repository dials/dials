/*
 * local.cc
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
#include <dials/algorithms/image/threshold/local.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  template <typename FloatType>
  void local_threshold_suite() {

    def("niblack", &niblack<FloatType>, (
      arg("image"),
      arg("size"),
      arg("n_sigma")));

    def("sauvola", &sauvola<FloatType>, (
      arg("image"),
      arg("size"),
      arg("k"),
      arg("r")));

    def("fano", &fano<FloatType>, (
      arg("image"),
      arg("size"),
      arg("n_sigma")));

    def("fano_masked", &fano_masked<FloatType>, (
      arg("image"),
      arg("mask"),
      arg("size"),
      arg("min_count"),
      arg("n_sigma")));

    def("gain", &gain<FloatType>, (
      arg("image"),
      arg("mask"),
      arg("gain"),
      arg("size"),
      arg("min_count"),
      arg("n_sigma")));

    def("kabsch", &kabsch<FloatType>, (
      arg("image"),
      arg("mask"),
      arg("size"),
      arg("n_sigma_b"),
      arg("n_sigma_s"),
      arg("min_count")));

    def("kabsch_w_gain", &kabsch_w_gain<FloatType>, (
      arg("image"),
      arg("mask"),
      arg("gain"),
      arg("size"),
      arg("n_sigma_b"),
      arg("n_sigma_s"),
      arg("min_count")));
  }

  void export_local() {
    local_threshold_suite<float>();
    local_threshold_suite<double>();


    class_<DispersionThreshold>("DispersionThreshold", no_init)
      .def(init< int2,
                 int2,
                 double,
                 double,
                 double,
                 int >())
      .def("__call__", &DispersionThreshold::threshold<int>)
      .def("__call__", &DispersionThreshold::threshold<double>)
      .def("__call__", &DispersionThreshold::threshold_w_gain<int>)
      .def("__call__", &DispersionThreshold::threshold_w_gain<double>)
      ;


    class_<KabschDebug>("KabschDebug", no_init)
      .def(init<const af::const_ref<double, af::c_grid<2> > &,
                const af::const_ref<bool, af::c_grid<2> > &,
                int2, double, double, double, int>((
                    arg("image"),
                    arg("mask"),
                    arg("size"),
                    arg("n_sigma_b"),
                    arg("n_sigma_s"),
                    arg("threshold"),
                    arg("min_count"))))
      .def(init<const af::const_ref<double, af::c_grid<2> > &,
                const af::const_ref<bool, af::c_grid<2> > &,
                const af::const_ref<double, af::c_grid<2> > &,
                int2, double, double, double, int>((
                    arg("image"),
                    arg("mask"),
                    arg("gain"),
                    arg("size"),
                    arg("n_sigma_b"),
                    arg("n_sigma_s"),
                    arg("threshold"),
                    arg("min_count"))))
      .def("mean", &KabschDebug::mean)
      .def("variance", &KabschDebug::variance)
      .def("coefficient_of_variation", &KabschDebug::coefficient_of_variation)
      .def("global_mask", &KabschDebug::global_mask)
      .def("cv_mask", &KabschDebug::cv_mask)
      .def("value_mask", &KabschDebug::value_mask)
      .def("final_mask", &KabschDebug::final_mask)
      ;

  }

}}} // namespace = dials::algorithms::boost_python
