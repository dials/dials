/*
 * poisson_discriminator.cc
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
#include <dials/algorithms/background/poisson_discriminator.h>

namespace dials { namespace algorithms { namespace background {
  namespace boost_python {

  using namespace boost::python;

  template <typename FloatType>
  void poisson_discriminator_suite() {

    def("moment",
      &moment<FloatType>, (
        arg("data"),
        arg("centre"),
        arg("k")));

    def("is_poisson_distributed",
      &is_poisson_distributed<FloatType>, (
        arg("data"),
        arg("n_sigma")));
  }

  template <typename FloatType>
  struct PDOperator {
    typedef void (PoissonDiscriminator::*with_mask)(
      const af::const_ref<FloatType>&, af::ref<int>) const;
    typedef af::shared<int> (PoissonDiscriminator::*without_mask)(
      const af::const_ref<FloatType>&) const;
  };

  void export_poisson_discriminator()
  {
    poisson_discriminator_suite<float>();
    poisson_discriminator_suite<double>();

    PDOperator<float>::with_mask call_with_mask_float =
      &PoissonDiscriminator::operator()<float>;
    PDOperator<double>::with_mask call_with_mask_double =
      &PoissonDiscriminator::operator()<double>;
    PDOperator<float>::without_mask call_without_mask_float =
      &PoissonDiscriminator::operator()<float>;
    PDOperator<double>::without_mask call_without_mask_double =
      &PoissonDiscriminator::operator()<double>;

    class_<PoissonDiscriminator>(
        "PoissonDiscriminator")
      .def(init<std::size_t, double>((
        arg("min_data") = 10,
        arg("n_sigma") = 3.0)))
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
