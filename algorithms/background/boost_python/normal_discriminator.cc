/*
 * normal_discriminator.cc
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
#include <dials/algorithms/background/normal_discriminator.h>

namespace dials { namespace algorithms { namespace background {
  namespace boost_python {

  using namespace boost::python;

  template <typename FloatType>
  bool is_normally_distributed_wrapper(
      const af::const_ref<FloatType> &data,
      double n_sigma) {
    if (n_sigma <= 0) {
      return is_normally_distributed(data);
    }
    return is_normally_distributed(data, n_sigma);
  }

  template <typename FloatType>
  void normal_discriminator_suite() {

    def("maximum_n_sigma",
      &maximum_n_sigma<FloatType>, (
        arg("data")));

    def("minimum_n_sigma",
      &minimum_n_sigma<FloatType>, (
        arg("data")));

    def("absolute_maximum_n_sigma",
      &maximum_n_sigma<FloatType>, (
        arg("data")));

    def("is_normally_distributed",
      &is_normally_distributed_wrapper<FloatType>, (
        arg("data"),
        arg("n_sigma") = -1));
  }

  template <typename FloatType>
  struct NDOperator {
    typedef void (NormalDiscriminator::*with_mask)(
      const af::const_ref<FloatType>&, af::ref<int>) const;
    typedef af::shared<int> (NormalDiscriminator::*without_mask)(
      const af::const_ref<FloatType>&) const;
  };

  void export_normal_discriminator()
  {
    def("normal_expected_n_sigma",
      &normal_expected_n_sigma, (
        arg("n_obs")));

    normal_discriminator_suite<float>();
    normal_discriminator_suite<double>();

    NDOperator<float>::with_mask call_with_mask_float =
      &NormalDiscriminator::operator()<float>;
    NDOperator<double>::with_mask call_with_mask_double =
      &NormalDiscriminator::operator()<double>;
    NDOperator<float>::without_mask call_without_mask_float =
      &NormalDiscriminator::operator()<float>;
    NDOperator<double>::without_mask call_without_mask_double =
      &NormalDiscriminator::operator()<double>;

    class_<NormalDiscriminator>("NormalDiscriminator", no_init)
      .def(init<std::size_t>((
        arg("min_data") = 10)))
      .def("__call__",
        call_with_mask_float, (
          arg("shoebox"),
          arg("mask")))
      .def("__call__",
        call_with_mask_double, (
          arg("shoebox"),
          arg("mask")))
      .def("__call__",
        call_without_mask_float, (
          arg("shoebox")))
      .def("__call__",
        call_without_mask_double, (
          arg("shoebox")));
  }

}}}} // namespace = dials::algorithms::background::boost_python
