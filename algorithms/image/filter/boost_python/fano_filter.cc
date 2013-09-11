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

  void export_fano_filter() 
  {
    class_<FanoFilter>("FanoFilter", no_init)
      .def(init<const af::const_ref< double, af::c_grid<2> >&, int2>((
        arg("image"), 
        arg("size"))))
      .def("fano", &FanoFilter::fano)
      .def("mean", &FanoFilter::mean)
      .def("sample_variance", &FanoFilter::sample_variance)
      .def("mask", &FanoFilter::mask);

    class_<FanoFilterMasked>("FanoFilterMasked", no_init)
      .def(init<const af::const_ref< double, af::c_grid<2> >&, 
                const af::const_ref< int, af::c_grid<2> >&, int2, int>((
        arg("image"), 
        arg("mask"), 
        arg("size"), 
        arg("min_count"))))
      .def("fano", &FanoFilterMasked::fano)
      .def("mean", &FanoFilterMasked::mean)
      .def("sample_variance", &FanoFilterMasked::sample_variance)
      .def("mask", &FanoFilterMasked::mask)
      .def("count", &FanoFilterMasked::count);
  }

}}} // namespace = dials::algorithms::boost_python
