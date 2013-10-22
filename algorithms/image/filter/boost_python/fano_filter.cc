/*
 * fano_filter.cc
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
#include <dials/algorithms/image/filter/fano_filter.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  template <typename FloatType>
  void fano_filter_wrapper(const char *name) {
  
    typedef FanoFilter<FloatType> FanoFilterType;
  
    class_<FanoFilterType>(name, no_init)
      .def(init<const af::const_ref< FloatType, af::c_grid<2> >&, int2>((
        arg("image"), 
        arg("size"))))
      .def("fano", &FanoFilterType::fano)
      .def("mean", &FanoFilterType::mean)
      .def("sample_variance", &FanoFilterType::sample_variance);
  }

  template <typename FloatType>
  void fano_filter_masked_wrapper(const char *name) {
  
    typedef FanoFilterMasked<FloatType> FanoFilterType;
  
    class_<FanoFilterType>(name, no_init)
      .def(init<const af::const_ref< FloatType, af::c_grid<2> >&, 
                const af::const_ref< int, af::c_grid<2> >&, int2, int>((
        arg("image"), 
        arg("mask"), 
        arg("size"), 
        arg("min_count"))))
      .def("fano", &FanoFilterType::fano)
      .def("mean", &FanoFilterType::mean)
      .def("sample_variance", &FanoFilterType::sample_variance)
      .def("mask", &FanoFilterType::mask)
      .def("count", &FanoFilterType::count);
  }
  
  template <typename FloatType>
  FanoFilter<FloatType> make_fano_filter(
    const af::const_ref<FloatType, af::c_grid<2> > &image, int2 size) {
    return FanoFilter<FloatType>(image, size);
  }

  template <typename FloatType>
  FanoFilterMasked<FloatType> make_fano_filter_masked(
    const af::const_ref<FloatType, af::c_grid<2> > &image, 
    const af::const_ref<int, af::c_grid<2> > &mask, 
    int2 size, int min_size) {
    return FanoFilterMasked<FloatType>(image, mask, size, min_size);
  }
  
  template <typename FloatType>
  void fano_filter_suite() {
    def("fano_filter", &make_fano_filter<FloatType>, (
      arg("image"),
      arg("kernel")));
      
    def("fano_filter", &make_fano_filter_masked<FloatType>, (
      arg("image"),
      arg("mask"),
      arg("kernel"),
      arg("min_count")));  
  }

  void export_fano_filter() {
    fano_filter_wrapper<float>("FanoFilterFloat");
    fano_filter_wrapper<double>("FanoFilterDouble");
    fano_filter_masked_wrapper<float>("FanoFilterMaskedFloat");
    fano_filter_masked_wrapper<double>("FanoFilterMaskedDouble");
    
    fano_filter_suite<float>();
    fano_filter_suite<double>();
  }

}}} // namespace = dials::algorithms::boost_python
