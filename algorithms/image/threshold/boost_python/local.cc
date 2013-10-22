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

  template <typename FloatType>
  void local_threshold_suite() {
  
    def("niblack", &niblack<FloatType>, (
      arg("image"),
      arg("size"),
      arg("n_sigma")));
  
    def("sauvola", &sauvola<FloatType>, (
      arg("image"),
      arg("size"),
      arg("k"),
      arg("r")));
  
    def("fano", &fano<FloatType>, (
      arg("image"),
      arg("size"),
      arg("n_sigma")));
      
    def("fano_masked", &fano_masked<FloatType>, (
      arg("image"),
      arg("mask"),
      arg("size"),
      arg("min_count"),
      arg("n_sigma")));

    def("gain", &gain<FloatType>, (
      arg("image"),
      arg("mask"),
      arg("gain"),
      arg("size"),
      arg("min_count"),
      arg("n_sigma")));

    def("kabsch", &kabsch<FloatType>, (
      arg("image"),
      arg("mask"),
      arg("size"),
      arg("n_sigma_b"),
      arg("n_sigma_s")));
      
    def("kabsch_w_gain", &kabsch_w_gain<FloatType>, (
      arg("image"),
      arg("mask"),
      arg("gain"),
      arg("size"),
      arg("n_sigma_b"),
      arg("n_sigma_s")));  
  }

  void export_local() {
    local_threshold_suite<float>();
    local_threshold_suite<double>();
  }

}}} // namespace = dials::algorithms::boost_python
