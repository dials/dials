/*
 * index_of_dispersion_filter.cc
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
#include <dials/algorithms/image/filter/index_of_dispersion_filter.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  template <typename FloatType>
  void index_of_dispersion_filter_wrapper(const char *name) {
    typedef IndexOfDispersionFilter<FloatType> IndexOfDispersionFilterType;

    class_<IndexOfDispersionFilterType>(name, no_init)
      .def(init<const af::const_ref<FloatType, af::c_grid<2> > &, int2>(
        (arg("image"), arg("size"))))
      .def("index_of_dispersion", &IndexOfDispersionFilterType::index_of_dispersion)
      .def("mean", &IndexOfDispersionFilterType::mean)
      .def("sample_variance", &IndexOfDispersionFilterType::sample_variance);
  }

  template <typename FloatType>
  void index_of_dispersion_filter_masked_wrapper(const char *name) {
    typedef IndexOfDispersionFilterMasked<FloatType> IndexOfDispersionFilterType;

    class_<IndexOfDispersionFilterType>(name, no_init)
      .def(init<const af::const_ref<FloatType, af::c_grid<2> > &,
                const af::const_ref<int, af::c_grid<2> > &,
                int2,
                int>((arg("image"), arg("mask"), arg("size"), arg("min_count"))))
      .def("index_of_dispersion", &IndexOfDispersionFilterType::index_of_dispersion)
      .def("mean", &IndexOfDispersionFilterType::mean)
      .def("sample_variance", &IndexOfDispersionFilterType::sample_variance)
      .def("mask", &IndexOfDispersionFilterType::mask)
      .def("count", &IndexOfDispersionFilterType::count);
  }

  template <typename FloatType>
  IndexOfDispersionFilter<FloatType> make_index_of_dispersion_filter(
    const af::const_ref<FloatType, af::c_grid<2> > &image,
    int2 size) {
    return IndexOfDispersionFilter<FloatType>(image, size);
  }

  template <typename FloatType>
  IndexOfDispersionFilterMasked<FloatType> make_index_of_dispersion_filter_masked(
    const af::const_ref<FloatType, af::c_grid<2> > &image,
    const af::const_ref<int, af::c_grid<2> > &mask,
    int2 size,
    int min_size) {
    return IndexOfDispersionFilterMasked<FloatType>(image, mask, size, min_size);
  }

  template <typename FloatType>
  void index_of_dispersion_filter_suite() {
    def("index_of_dispersion_filter",
        &make_index_of_dispersion_filter<FloatType>,
        (arg("image"), arg("kernel")));

    def("index_of_dispersion_filter",
        &make_index_of_dispersion_filter_masked<FloatType>,
        (arg("image"), arg("mask"), arg("kernel"), arg("min_count")));
  }

  void export_index_of_dispersion_filter() {
    index_of_dispersion_filter_wrapper<float>("IndexOfDispersionFilterFloat");
    index_of_dispersion_filter_wrapper<double>("IndexOfDispersionFilterDouble");
    index_of_dispersion_filter_masked_wrapper<float>(
      "IndexOfDispersionFilterMaskedFloat");
    index_of_dispersion_filter_masked_wrapper<double>(
      "IndexOfDispersionFilterMaskedDouble");

    index_of_dispersion_filter_suite<float>();
    index_of_dispersion_filter_suite<double>();
  }

}}}  // namespace dials::algorithms::boost_python
