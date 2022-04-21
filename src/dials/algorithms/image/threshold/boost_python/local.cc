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
    def("niblack", &niblack<FloatType>, (arg("image"), arg("size"), arg("n_sigma")));

    def(
      "sauvola", &sauvola<FloatType>, (arg("image"), arg("size"), arg("k"), arg("r")));

    def("index_of_dispersion",
        &index_of_dispersion<FloatType>,
        (arg("image"), arg("size"), arg("n_sigma")));

    def("index_of_dispersion_masked",
        &index_of_dispersion_masked<FloatType>,
        (arg("image"), arg("mask"), arg("size"), arg("min_count"), arg("n_sigma")));

    def("gain",
        &gain<FloatType>,
        (arg("image"),
         arg("mask"),
         arg("gain"),
         arg("size"),
         arg("min_count"),
         arg("n_sigma")));

    def("dispersion",
        &dispersion<FloatType>,
        (arg("image"),
         arg("mask"),
         arg("size"),
         arg("n_sigma_b"),
         arg("n_sigma_s"),
         arg("min_count")));

    def("dispersion_w_gain",
        &dispersion_w_gain<FloatType>,
        (arg("image"),
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
      .def(init<int2, int2, double, double, double, int>())
      .def("__call__", &DispersionThreshold::threshold<int>)
      .def("__call__", &DispersionThreshold::threshold<double>)
      .def("__call__", &DispersionThreshold::threshold_w_gain<int>)
      .def("__call__", &DispersionThreshold::threshold_w_gain<double>);

    class_<DispersionThresholdDebug>("DispersionThresholdDebug", no_init)
      .def(init<const af::const_ref<double, af::c_grid<2> > &,
                const af::const_ref<bool, af::c_grid<2> > &,
                int2,
                double,
                double,
                double,
                int>((arg("image"),
                      arg("mask"),
                      arg("size"),
                      arg("n_sigma_b"),
                      arg("n_sigma_s"),
                      arg("threshold"),
                      arg("min_count"))))
      .def(init<const af::const_ref<double, af::c_grid<2> > &,
                const af::const_ref<bool, af::c_grid<2> > &,
                const af::const_ref<double, af::c_grid<2> > &,
                int2,
                double,
                double,
                double,
                int>((arg("image"),
                      arg("mask"),
                      arg("gain"),
                      arg("size"),
                      arg("n_sigma_b"),
                      arg("n_sigma_s"),
                      arg("threshold"),
                      arg("min_count"))))
      .def("mean", &DispersionThresholdDebug::mean)
      .def("variance", &DispersionThresholdDebug::variance)
      .def("index_of_dispersion", &DispersionThresholdDebug::index_of_dispersion)
      .def("global_mask", &DispersionThresholdDebug::global_mask)
      .def("cv_mask", &DispersionThresholdDebug::cv_mask)
      .def("value_mask", &DispersionThresholdDebug::value_mask)
      .def("final_mask", &DispersionThresholdDebug::final_mask);

    class_<DispersionExtendedThresholdDebug>("DispersionExtendedThresholdDebug",
                                             no_init)
      .def(init<const af::const_ref<double, af::c_grid<2> > &,
                const af::const_ref<bool, af::c_grid<2> > &,
                int2,
                double,
                double,
                double,
                int>((arg("image"),
                      arg("mask"),
                      arg("size"),
                      arg("n_sigma_b"),
                      arg("n_sigma_s"),
                      arg("threshold"),
                      arg("min_count"))))
      .def(init<const af::const_ref<double, af::c_grid<2> > &,
                const af::const_ref<bool, af::c_grid<2> > &,
                const af::const_ref<double, af::c_grid<2> > &,
                int2,
                double,
                double,
                double,
                int>((arg("image"),
                      arg("mask"),
                      arg("gain"),
                      arg("size"),
                      arg("n_sigma_b"),
                      arg("n_sigma_s"),
                      arg("threshold"),
                      arg("min_count"))))
      .def("mean", &DispersionExtendedThresholdDebug::mean)
      .def("variance", &DispersionExtendedThresholdDebug::variance)
      .def("index_of_dispersion",
           &DispersionExtendedThresholdDebug::index_of_dispersion)
      .def("global_mask", &DispersionExtendedThresholdDebug::global_mask)
      .def("cv_mask", &DispersionExtendedThresholdDebug::cv_mask)
      .def("value_mask", &DispersionExtendedThresholdDebug::value_mask)
      .def("final_mask", &DispersionExtendedThresholdDebug::final_mask);

    class_<DispersionExtendedThreshold>("DispersionExtendedThreshold", no_init)
      .def(init<int2, int2, double, double, double, int>())
      /* .def("__call__", &DispersionExtendedThreshold::threshold<int>) */
      .def("__call__", &DispersionExtendedThreshold::threshold<double>)
      /* .def("__call__", &DispersionExtendedThreshold::threshold_w_gain<int>) */
      .def("__call__", &DispersionExtendedThreshold::threshold_w_gain<double>);
  }

}}}  // namespace dials::algorithms::boost_python
