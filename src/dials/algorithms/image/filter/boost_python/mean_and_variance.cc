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

  template <typename FloatType>
  void mean_and_variance_filter_wrapper(const char *name) {
    typedef MeanAndVarianceFilter<FloatType> MeanAndVarianceFilterType;

    class_<MeanAndVarianceFilterType>(name, no_init)
      .def(init<const af::const_ref<FloatType, af::c_grid<2> > &, int2>(
        (arg("image"), arg("size"))))
      .def("mean", &MeanAndVarianceFilterType::mean)
      .def("variance", &MeanAndVarianceFilterType::variance)
      .def("sample_variance", &MeanAndVarianceFilterType::sample_variance);
  }

  template <typename FloatType>
  void mean_and_variance_filter_masked_wrapper(const char *name) {
    typedef MeanAndVarianceFilterMasked<FloatType> MeanAndVarianceFilterType;

    class_<MeanAndVarianceFilterType>(name, no_init)
      .def(init<const af::const_ref<FloatType, af::c_grid<2> > &,
                const af::const_ref<int, af::c_grid<2> > &,
                int2,
                int>((arg("image"), arg("mask"), arg("size"), arg("min_size"))))
      .def("mean", &MeanAndVarianceFilterType::mean)
      .def("variance", &MeanAndVarianceFilterType::variance)
      .def("sample_variance", &MeanAndVarianceFilterType::sample_variance)
      .def("mask", &MeanAndVarianceFilterType::mask)
      .def("count", &MeanAndVarianceFilterType::count);
  }

  template <typename FloatType>
  MeanAndVarianceFilter<FloatType> make_mean_and_variance_filter(
    const af::const_ref<FloatType, af::c_grid<2> > &image,
    int2 size) {
    return MeanAndVarianceFilter<FloatType>(image, size);
  }

  template <typename FloatType>
  MeanAndVarianceFilterMasked<FloatType> make_mean_and_variance_filter_masked(
    const af::const_ref<FloatType, af::c_grid<2> > &image,
    const af::const_ref<int, af::c_grid<2> > &mask,
    int2 size,
    int min_size) {
    return MeanAndVarianceFilterMasked<FloatType>(image, mask, size, min_size);
  }

  template <typename FloatType>
  void mean_and_variance_filter_suite() {
    def("mean_filter", &mean_filter<FloatType>, (arg("image"), arg("size")));

    def("mean_filter",
        &mean_filter_masked<FloatType>,
        (arg("image"),
         arg("mask"),
         arg("size"),
         arg("min_count"),
         arg("ignore_masked") = true));

    def("mean_and_variance_filter",
        &make_mean_and_variance_filter<FloatType>,
        (arg("image"), arg("kernel")));

    def("mean_and_variance_filter",
        &make_mean_and_variance_filter_masked<FloatType>,
        (arg("image"), arg("mask"), arg("kernel"), arg("min_count")));
  }

  void export_mean_and_variance() {
    mean_and_variance_filter_wrapper<float>("MeanAndVarianceFilterFloat");
    mean_and_variance_filter_wrapper<double>("MeanAndVarianceFilterDouble");
    mean_and_variance_filter_masked_wrapper<float>("MeanAndVarianceFilterMaskedFloat");
    mean_and_variance_filter_masked_wrapper<double>(
      "MeanAndVarianceFilterMaskedDouble");

    mean_and_variance_filter_suite<float>();
    mean_and_variance_filter_suite<double>();
  }

}}}  // namespace dials::algorithms::boost_python
