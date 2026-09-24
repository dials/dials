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
#include <boost/python/numpy.hpp>
#include <boost/python/def.hpp>
#include <dials/algorithms/background/modeller.h>

namespace np = boost::python::numpy;

namespace dials { namespace algorithms { namespace background { namespace boost_python {

  using namespace boost::python;

  struct BSPickleSuite : pickle_suite {
    static tuple getinitargs(BackgroundStatistics& obj) {
      return make_tuple(ImageVolume<>(0, 1, 2, 3));
    }

    static object getstate(BackgroundStatistics& obj) {
      return obj.get_all();
    }

    static void setstate(BackgroundStatistics& obj, object state) {
      obj.set_all(state);
    }
  };

  struct MPBSPickleSuite : pickle_suite {
    static tuple getinitargs(MultiPanelBackgroundStatistics& obj) {
      return make_tuple(MultiPanelImageVolume<>());
    }

    static object getstate(MultiPanelBackgroundStatistics& obj) {
      return make_tuple(obj.get_statistics());
    }

    static void setstate(MultiPanelBackgroundStatistics& obj, object state) {
      obj.set_statistics(
        extract<af_shared_with_getitem<BackgroundStatistics>>(state[0]));
    }
  };

  struct SBSPickleSuite : pickle_suite {
    static tuple getinitargs(af_shared_with_getitem<BackgroundStatistics>& obj) {
      return tuple();
    }

    static object getstate(af_shared_with_getitem<BackgroundStatistics>& obj) {
      list arg_lst;
      for (int i = 0; i < obj.size(); i++) {
        arg_lst.append(obj[i]);
      }
      return tuple(arg_lst);
    }

    static void setstate(af_shared_with_getitem<BackgroundStatistics>& obj,
                         object state) {
      for (int i = 0; i < len(state); i++) {
        obj.push_back(extract<BackgroundStatistics>(state[i]));
      }
    }
  };

  BOOST_PYTHON_MODULE(dials_algorithms_background_modeller_ext) {
    np::initialize();

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
      .def("mask", &BackgroundStatistics::mask)
      .def_pickle(BSPickleSuite());

    class_<MultiPanelBackgroundStatistics>("MultiPanelBackgroundStatistics", no_init)
      .def(init<const MultiPanelImageVolume<>&>())
      .def("get", &MultiPanelBackgroundStatistics::get)
      .def("__len__", &MultiPanelBackgroundStatistics::size)
      .def("__iadd__", &MultiPanelBackgroundStatistics::operator+=)
      .def("get_statistics", &MultiPanelBackgroundStatistics::get_statistics)
      .def_pickle(MPBSPickleSuite());

    class_<af_shared_with_getitem<BackgroundStatistics>>("shared_BackgroundStatistics")
      .def_pickle(SBSPickleSuite());
  }

}}}}  // namespace dials::algorithms::background::boost_python
