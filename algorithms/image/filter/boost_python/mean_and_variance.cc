/*
 * mean_and_variance.cc
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
#include <dials/algorithms/image/filter/mean_and_variance.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_mean_and_variance() 
  {
    def("mean_filter", &mean_filter, (
      arg("image"), 
      arg("size")));
    def("mean_filter_masked", &mean_filter_masked, (
      arg("image"), 
      arg("mask"), 
      arg("size"), 
      arg("min_count")));

    class_<MeanAndVarianceFilter>("MeanAndVarianceFilter", no_init)
      .def(init<const af::const_ref<double, af::c_grid<2> >&, int2>((
          arg("image"), 
          arg("size"))))
      .def("mean", &MeanAndVarianceFilter::mean)
      .def("variance", &MeanAndVarianceFilter::variance)
      .def("sample_variance", &MeanAndVarianceFilter::sample_variance);

    class_<MeanAndVarianceFilterMasked>("MeanAndVarianceFilterMasked", no_init)
      .def(init<const af::const_ref<double, af::c_grid<2> >&, 
                const af::const_ref<int, af::c_grid<2> >&, int2, int>((
          arg("image"), 
          arg("mask"), 
          arg("size"), 
          arg("min_size"))))
      .def("mean", &MeanAndVarianceFilterMasked::mean)
      .def("variance", &MeanAndVarianceFilterMasked::variance)
      .def("sample_variance", &MeanAndVarianceFilterMasked::sample_variance)
      .def("mask", &MeanAndVarianceFilterMasked::mask)
      .def("count", &MeanAndVarianceFilterMasked::count);
  }

}}} // namespace = dials::algorithms::boost_python
