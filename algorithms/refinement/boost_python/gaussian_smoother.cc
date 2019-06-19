/*
 * gaussian_smoother.cc
 *
 *  Copyright (C) (2016) STFC Rutherford Appleton Laboratory, UK.
 *
 *  Author: David Waterman.
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../gaussian_smoother.h"

using namespace boost::python;

namespace dials { namespace refinement { namespace boost_python {

  void export_gaussian_smoother() {
    class_<GaussianSmoother>("GaussianSmoother", no_init)
      .def(init<vec2<double>, std::size_t>((arg("x_range"), arg("num_intervals"))))
      .def("set_smoothing", &GaussianSmoother::set_smoothing)
      .def("num_values", &GaussianSmoother::num_values)
      .def("num_samples", &GaussianSmoother::num_samples)
      .def("num_average", &GaussianSmoother::num_average)
      .def("sigma", &GaussianSmoother::sigma)
      .def("spacing", &GaussianSmoother::spacing)
      .def("positions", &GaussianSmoother::positions)
      .def("value_weight", &GaussianSmoother::value_weight)
      .def("multi_value_weight", &GaussianSmoother::multi_value_weight);

    class_<SingleValueWeights>("SingleValueWeights", no_init)
      .def("get_value", &SingleValueWeights::get_value)
      .def("get_weight", &SingleValueWeights::get_weight)
      .def("get_sumweight", &SingleValueWeights::get_sumweight);

    class_<MultiValueWeights>("MultiValueWeights", no_init)
      .def("get_value", &MultiValueWeights::get_value)
      .def("get_weight", &MultiValueWeights::get_weight)
      .def("get_sumweight", &MultiValueWeights::get_sumweight);
  }

}}}  // namespace dials::refinement::boost_python
