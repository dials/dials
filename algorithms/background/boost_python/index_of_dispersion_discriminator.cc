/*
 * index_of_dispersion_discriminator.cc
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
#include <dials/algorithms/background/index_of_dispersion_discriminator.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_index_of_dispersion_discriminator()
  {

    // Export normality test
    def("is_poisson_distributed_iod", &is_poisson_distributed_iod, (
      arg("data"), arg("n_sigma") = 3));  
    
    // Overloads for call method
    void (IndexOfDispersionDiscriminator::*call_shoebox_and_mask)(
      const flex_int&, flex_int &) const = 
        &IndexOfDispersionDiscriminator::operator();
    flex_int (IndexOfDispersionDiscriminator::*call_shoebox)(const flex_int&) 
      const = &IndexOfDispersionDiscriminator::operator();
  
    class_<IndexOfDispersionDiscriminator, bases<DiscriminatorStrategy> >(
        "IndexOfDispersionDiscriminator")
      .def(init<std::size_t, double>((
        arg("min_data") = 10,
        arg("n_sigma") = 3.0)))
      .def("__call__", call_shoebox_and_mask, (
        arg("shoebox"),
        arg("mask")))
      .def("__call__", call_shoebox, (
        arg("shoebox")));
  }

}}} // namespace = dials::algorithms::boost_python
