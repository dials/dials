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
#include <dials/algorithms/background/outlier_rejector.h>
#include <dials/algorithms/background/normal_outlier_rejector.h>
#include <dials/algorithms/background/truncated_outlier_rejector.h>
#include <dials/algorithms/background/nsigma_outlier_rejector.h>

namespace dials { namespace algorithms { namespace background {
  namespace boost_python {

  using namespace boost::python;

  void export_outlier_rejector()
  {
    class_<OutlierRejector, boost::noncopyable>("OutlierRejector", no_init)
      .def("__call__", &OutlierRejector::mark, (
            arg("data"),
            arg("mask")));

    class_<TruncatedOutlierRejector, bases<OutlierRejector> >(
        "TruncatedOutlierRejector", no_init)
      .def(init<double, double>((
        arg("lower") = 0.01,
        arg("upper") = 0.01)));

    class_<NSigmaOutlierRejector, bases<OutlierRejector> >(
        "NSigmaOutlierRejector", no_init)
      .def(init<double, double>((
        arg("lower") = 3.0,
        arg("upper") = 3.0)));

    class_<NormalOutlierRejector, bases<OutlierRejector> >(
        "NormalOutlierRejector", no_init)
      .def(init<std::size_t>((
        arg("min_data") = 10)));
  }

}}}} // namespace = dials::algorithms::background::boost_python
