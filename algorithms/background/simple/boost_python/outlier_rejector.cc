/*
 * outlier_rejector.cc
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
#include <dials/algorithms/background/simple/outlier_rejector.h>
#include <dials/algorithms/background/simple/normal_outlier_rejector.h>
#include <dials/algorithms/background/simple/truncated_outlier_rejector.h>
#include <dials/algorithms/background/simple/nsigma_outlier_rejector.h>
#include <dials/algorithms/background/simple/mosflm_outlier_rejector.h>
#include <dials/algorithms/background/simple/tukey_outlier_rejector.h>

namespace dials { namespace algorithms { namespace background { namespace boost_python {

  using namespace boost::python;

  template <typename FloatType>
  bool is_normally_distributed_wrapper(const af::const_ref<FloatType> &data,
                                       double n_sigma) {
    if (n_sigma <= 0) {
      return is_normally_distributed(data);
    }
    return is_normally_distributed(data, n_sigma);
  }

  template <typename FloatType>
  void normal_discriminator_suite() {
    def("maximum_n_sigma", &maximum_n_sigma<FloatType>, (boost::python::arg("data")));

    def("minimum_n_sigma", &minimum_n_sigma<FloatType>, (boost::python::arg("data")));

    def("absolute_maximum_n_sigma",
        &maximum_n_sigma<FloatType>,
        (boost::python::arg("data")));

    def("is_normally_distributed",
        &is_normally_distributed_wrapper<FloatType>,
        (boost::python::arg("data"), boost::python::arg("n_sigma") = -1));
  }

  void export_outlier_rejector() {
    def("normal_expected_n_sigma",
        &normal_expected_n_sigma,
        (boost::python::arg("n_obs")));

    normal_discriminator_suite<float>();
    normal_discriminator_suite<double>();

    class_<OutlierRejector, boost::noncopyable>("OutlierRejector", no_init)
      .def("__call__",
           &OutlierRejector::mark,
           (boost::python::arg("data"), boost::python::arg("mask")));

    class_<TruncatedOutlierRejector, bases<OutlierRejector> >(
      "TruncatedOutlierRejector", no_init)
      .def(init<double, double>(
        (boost::python::arg("lower") = 0.01, boost::python::arg("upper") = 0.01)));

    class_<NSigmaOutlierRejector, bases<OutlierRejector> >("NSigmaOutlierRejector",
                                                           no_init)
      .def(init<double, double>(
        (boost::python::arg("lower") = 3.0, boost::python::arg("upper") = 3.0)));

    class_<NormalOutlierRejector, bases<OutlierRejector> >("NormalOutlierRejector",
                                                           no_init)
      .def(init<std::size_t>((boost::python::arg("min_data") = 10)));

    class_<MosflmOutlierRejector, bases<OutlierRejector> >("MosflmOutlierRejector",
                                                           no_init)
      .def(init<double, double>(
        (boost::python::arg("fraction") = 1.0, boost::python::arg("n_sigma") = 3)));

    class_<TukeyOutlierRejector, bases<OutlierRejector> >("TukeyOutlierRejector",
                                                          no_init)
      .def(init<double, double>(
        (boost::python::arg("lower") = 1.5, boost::python::arg("upper") = 1.5)));
  }

}}}}  // namespace dials::algorithms::background::boost_python
