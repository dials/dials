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

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  bool is_normally_distributed_wrapper(const const_ref<double> &data, 
      double n_sigma) {
    if (n_sigma <= 0) {
      return is_normally_distributed(data);
    }
    return is_normally_distributed(data, n_sigma);
  }

  void export_normal_discriminator()
  {
    // Export the expected number of sdevs
    def("normal_expected_n_sigma", &normal_expected_n_sigma, (arg("n_obs")));
  
    // Export maximum_n_sigma
    def("maximum_n_sigma", &maximum_n_sigma, (arg("data")));
  
    // Export minimum_n_sigma
    def("minimum_n_sigma", &minimum_n_sigma, (arg("data")));
  
    // Export maximum_n_sigma
    def("absolute_maximum_n_sigma", &maximum_n_sigma, (arg("data")));
  
    // Export normality test
    def("is_normally_distributed", &is_normally_distributed_wrapper, (
      arg("data"), arg("n_sigma") = -1));  
    
    // Overloads for call method
    void (NormalDiscriminator::*call_shoebox_and_mask)(const flex_double&,
        flex_int &) const = &NormalDiscriminator::operator();
    flex_int (NormalDiscriminator::*call_shoebox)(const flex_double&) const =
      &NormalDiscriminator::operator();
    void (NormalDiscriminator::*call_reflection)(Reflection &) const =
      &NormalDiscriminator::operator();  
  
    class_<NormalDiscriminator>(
        "NormalDiscriminator")
      .def(init<std::size_t, double>((
        arg("min_data") = 10,
        arg("n_sigma") = 3.0)))
      .def("__call__", call_shoebox_and_mask, (
        arg("shoebox"),
        arg("mask")))
      .def("__call__", call_shoebox, (
        arg("shoebox")))
      .def("__call__", call_reflection, (
        arg("reflection")));
  }

}}} // namespace = dials::algorithms::boost_python
