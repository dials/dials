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

  struct BackgroundStatisticsPickleSuite : pickle_suite {
    static tuple getinitargs(const BackgroundStatistics& obj) {
      return make_tuple(obj.sum(), obj.sum_sq(), obj.num(), obj.min(), obj.max());
    }
  };

  struct MultiPanelBackgroundStatisticsPickleSuite : pickle_suite {
    static tuple getinitargs(const MultiPanelBackgroundStatistics& obj) {
      return tuple();
    }

    static tuple getstate(const MultiPanelBackgroundStatistics& obj) {
      list statistics;
      for (std::size_t i = 0; i < obj.size(); ++i) {
        statistics.append(obj.get(i));
      }
      return tuple(statistics);
    }

    static void setstate(MultiPanelBackgroundStatistics& obj, tuple state) {
      DIALS_ASSERT(obj.size() == 0);
      for (std::size_t i = 0; i < len(state); ++i) {
        obj.add(extract<const BackgroundStatistics&>(state[i]));
      }
    }
  };

  BOOST_PYTHON_MODULE(dials_algorithms_background_modeller_ext) {
    class_<BackgroundStatistics>("BackgroundStatistics", no_init)
      .def(init<const ImageVolume<>&>())
      .def(init<const af::const_ref<double, af::c_grid<2> >&,
                const af::const_ref<double, af::c_grid<2> >&,
                const af::const_ref<int, af::c_grid<2> >&,
                const af::const_ref<double, af::c_grid<2> >&,
                const af::const_ref<double, af::c_grid<2> >&>(
        (arg("sum"), arg("sum_sq"), arg("num"), arg("min"), arg("max"))))
      .def("sum", &BackgroundStatistics::sum)
      .def("sum_sq", &BackgroundStatistics::sum_sq)
      .def("num", &BackgroundStatistics::num)
      .def("min", &BackgroundStatistics::min)
      .def("max", &BackgroundStatistics::max)
      .def("__iadd__", &BackgroundStatistics::operator+=)
      .def("mean", &BackgroundStatistics::mean)
      .def("variance", &BackgroundStatistics::variance)
      .def("dispersion", &BackgroundStatistics::dispersion)
      .def("mask", &BackgroundStatistics::mask)
      .def_pickle(BackgroundStatisticsPickleSuite());

    class_<MultiPanelBackgroundStatistics>("MultiPanelBackgroundStatistics", no_init)
      .def(init<const MultiPanelImageVolume<>&>())
      .def(init<>())
      .def("get", &MultiPanelBackgroundStatistics::get)
      .def("__len__", &MultiPanelBackgroundStatistics::size)
      .def("__iadd__", &MultiPanelBackgroundStatistics::operator+=)
      .def_pickle(MultiPanelBackgroundStatisticsPickleSuite());
  }

}}}}  // namespace dials::algorithms::background::boost_python
