/*
 * reciprocal_space_transform.cc
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
#include <dials/algorithms/integration/reciprocal_space_transform.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_reciprocal_space_transform() 
  {
    class_<ReciprocalSpaceTransformE3Fraction>(
        "ReciprocalSpaceTransformE3Fraction", no_init)
      .def(init<const Scan&, double, double, int>((
        arg("scan"),
        arg("mosaicity"),
        arg("n_sigma"),
        arg("grid_size_e3"))))
      .def("__call__",
        &ReciprocalSpaceTransformE3Fraction::operator(), (
          arg("bbox_z"),
          arg("phi"),
          arg("zeta")));
  
    class_<ReciprocalSpaceTransformDetectorLabCoords>(
        "ReciprocalSpaceTransformDetectorLabCoords")
      .def("__call__",
        &ReciprocalSpaceTransformDetectorLabCoords::operator(), (
          arg("detector"),
          arg("scan"),
          arg("n_div")));
  
  
    flex_double (ReciprocalSpaceTransform::*call_single)(
      const flex_double&, const flex_int&, int6, vec3<double>, double) const = 
        &ReciprocalSpaceTransform::operator();
  
    void (ReciprocalSpaceTransform::*call_reflection)(
      Reflection&) const = &ReciprocalSpaceTransform::operator();
    
    void (ReciprocalSpaceTransform::*call_reflection_list)(
      ReflectionList&) const = &ReciprocalSpaceTransform::operator();
  
    class_<ReciprocalSpaceTransform>("ReciprocalSpaceTransform", no_init)
      .def(init<const Beam&,
                const Detector&,
                const Goniometer&,
                const Scan&,
                double,
                double,
                std::size_t,
                std::size_t>((
        arg("beam"),
        arg("detector"),
        arg("goniometer"),
        arg("scan"),
        arg("mosaicity"),
        arg("n_sigma"),
        arg("grid_size"),
        arg("n_div"))))
      .def("__call__", call_single, (
        arg("pixels"),
        arg("mask"),
        arg("bbox"),
        arg("s1"),
        arg("phi")))
      .def("__call__", call_reflection, (
        arg("reflection")))
      .def("__call__", call_reflection_list, (
        arg("reflection_list")));
  }

}}} // namespace = dials::integration::boost_python
