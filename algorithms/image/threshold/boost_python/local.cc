/*
 * local.cc
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
#include <dials/algorithms/image/threshold/local.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_local() 
  {
    def("niblack", &niblack, (
      arg("image"),
      arg("size"),
      arg("n_sigma")));
  
    def("sauvola", &sauvola, (
      arg("image"),
      arg("size"),
      arg("k"),
      arg("r")));
  
    def("fano", &fano, (
      arg("image"),
      arg("size"),
      arg("n_sigma")));
      
    def("fano_masked", &fano_masked, (
      arg("image"),
      arg("mask"),
      arg("size"),
      arg("min_count"),
      arg("n_sigma")));


    def("gain", &gain, (
      arg("image"),
      arg("mask"),
      arg("gain"),
      arg("size"),
      arg("min_count"),
      arg("n_sigma")));
    
  }

}}} // namespace = dials::algorithms::boost_python
