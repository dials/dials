/*
 * filter_ext.cc
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
#include <dials/algorithms/integration/filter.h>

namespace dials { namespace algorithms { namespace filter { 
    namespace boost_python {

  using namespace boost::python;

  void export_is_zeta_valid()
  {
    def("is_zeta_valid", 
      (bool(*)(vec3<double>, vec3<double>, vec3<double>, double))
      &is_zeta_valid, (
        arg("m2"), arg("s0"), arg("s1"), arg("zeta_min")));    
    def("is_zeta_valid",
      (bool(*)(const CoordinateSystem&, double)) 
      &is_zeta_valid, (
        arg("cs"), arg("zeta_min")));
    def("is_zeta_valid",
      (bool(*)(vec3<double>, vec3<double>, const Reflection&, double))
      &is_zeta_valid, (
        arg("m2"), arg("s0"), arg("r"), arg("zeta_min")));
    def("is_zeta_valid",
      (bool(*)(const Goniometer&, const Beam&, const Reflection&, double))
      &is_zeta_valid, (
        arg("g"), arg("b"), arg("r"), arg("zeta_min")));
  }

  void export_is_xds_small_angle_valid()
  {
    def("is_xds_small_angle_valid", 
      (bool(*)(vec3<double>, vec3<double>, vec3<double>, double))
      &is_xds_small_angle_valid, (
        arg("m2"), arg("s0"), arg("s1"), arg("delta_m")));    
    def("is_xds_small_angle_valid",
      (bool(*)(const CoordinateSystem&, double)) 
      &is_xds_small_angle_valid, (
        arg("cs"), arg("delta_m")));
    def("is_xds_small_angle_valid",
      (bool(*)(vec3<double>, vec3<double>, const Reflection&, double))
      &is_xds_small_angle_valid, (
        arg("m2"), arg("s0"), arg("r"), arg("delta_m")));
    def("is_xds_small_angle_valid",
      (bool(*)(const Goniometer&, const Beam&, const Reflection&, double))
      &is_xds_small_angle_valid, (
        arg("g"), arg("b"), arg("r"), arg("delta_m")));
  }
  
  void export_is_xds_angle_valid()
  {
    def("is_xds_angle_valid", 
      (bool(*)(vec3<double>, vec3<double>, vec3<double>, double))
      &is_xds_angle_valid, (
        arg("m2"), arg("s0"), arg("s1"), arg("delta_m")));    
    def("is_xds_angle_valid",
      (bool(*)(const CoordinateSystem&, double)) 
      &is_xds_angle_valid, (
        arg("cs"), arg("delta_m")));
    def("is_xds_angle_valid",
      (bool(*)(vec3<double>, vec3<double>, const Reflection&, double))
      &is_xds_angle_valid, (
        arg("m2"), arg("s0"), arg("r"), arg("delta_m")));
    def("is_xds_angle_valid",
      (bool(*)(const Goniometer&, const Beam&, const Reflection&, double))
      &is_xds_angle_valid, (
        arg("g"), arg("b"), arg("r"), arg("delta_m")));
  }

  BOOST_PYTHON_MODULE(dials_algorithms_integration_filter_ext)
  {
    export_is_zeta_valid();
    export_is_xds_small_angle_valid();
    export_is_xds_angle_valid();
  }

}}}} // namespace = dials::algorithms::filter::boost_python
