/*
 * background_intensity.cc
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
#include <dials/algorithms/integration/background_intensity.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  bool is_normally_distributed_wrapper(const const_ref<double> &data, 
      double n_sigma) {
    if (n_sigma <= 0) {
      return is_normally_distributed(data);
    }
    return is_normally_distributed(data, n_sigma);
  }

  void export_background_intensity() 
  {
    // Export the expected number of sdevs
    def("expected_n_sigma", &expected_n_sigma, (arg("n_obs")));
  
    // Export maximum_n_sigma
    def("maximum_n_sigma", &maximum_n_sigma, (arg("data")));
  
    // Export normality test
    def("is_normally_distributed", &is_normally_distributed_wrapper, (
      arg("data"), arg("n_sigma") = -1));
      
    // Export background intensity calculation
    def("background_intensity", &background_intensity, (
      arg("data"), arg("min_data") = 10, arg("n_sigma") = -1));
  }

}}} // namespace = dials::algorithms::boost_python
