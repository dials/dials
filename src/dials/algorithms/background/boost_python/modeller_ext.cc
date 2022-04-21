/*
 * modeller_ext.cc
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
#include <dials/algorithms/background/modeller.h>

namespace dials { namespace algorithms { namespace background { namespace boost_python {

  using namespace boost::python;

  BOOST_PYTHON_MODULE(dials_algorithms_background_modeller_ext) {
    class_<BackgroundStatistics>("BackgroundStatistics", no_init)
      .def(init<const ImageVolume<>&>())
      .def("sum", &BackgroundStatistics::sum)
      .def("sum_sq", &BackgroundStatistics::sum_sq)
      .def("num", &BackgroundStatistics::num)
      .def("min", &BackgroundStatistics::min)
      .def("max", &BackgroundStatistics::max)
      .def("__iadd__", &BackgroundStatistics::operator+=)
      .def("mean", &BackgroundStatistics::mean)
      .def("variance", &BackgroundStatistics::variance)
      .def("dispersion", &BackgroundStatistics::dispersion)
      .def("mask", &BackgroundStatistics::mask);

    class_<MultiPanelBackgroundStatistics>("MultiPanelBackgroundStatistics", no_init)
      .def(init<const MultiPanelImageVolume<>&>())
      .def("get", &MultiPanelBackgroundStatistics::get)
      .def("__len__", &MultiPanelBackgroundStatistics::size)
      .def("__iadd__", &MultiPanelBackgroundStatistics::operator+=);
  }

}}}}  // namespace dials::algorithms::background::boost_python
