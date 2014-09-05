/*
 * bbox_calculator.cc
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
#include <dials/algorithms/shoebox/bbox_calculator.h>

namespace dials { namespace algorithms { namespace shoebox {
  namespace boost_python {

  using namespace boost::python;

  void export_bbox_calculator()
  {
    int6 (BBoxCalculator::*calculate_single)(
      vec3 <double>, double, std::size_t) const =
        &BBoxCalculator::operator();
    af::shared<int6> (BBoxCalculator::*calculate_array) (
      const af::const_ref< vec3<double> >&,
      const af::const_ref<double> &,
      const af::const_ref<std::size_t>&) const =
        &BBoxCalculator::operator();

    class_ <BBoxCalculator> ("BBoxCalculator", no_init)
      .def(init <const Beam&,
                 const Detector&,
                 const Goniometer&,
                 const Scan&,
                 const af::const_ref<double>&,
                 const af::const_ref<double>&> ((
        arg("beam"),
        arg("detector"),
        arg("goniometer"),
        arg("scan"),
        arg("delta_divergence"),
        arg("delta_mosaicity"))))
      .def(init <const Beam&,
                 const Detector&,
                 const Goniometer&,
                 const Scan&,
                 double,
                 double > ((
        arg("beam"),
        arg("detector"),
        arg("goniometer"),
        arg("scan"),
        arg("delta_divergence"),
        arg("delta_mosaicity"))))
      .def("__call__", calculate_single, (
        arg("s1"), arg("phi"), arg("panel")))
      .def("__call__", calculate_array, (
        arg("s1"), arg("phi"), arg("panel")));

    class_ <BBoxMultiCalculator>("BBoxMultiCalculator")
      .def("append", &BBoxMultiCalculator::push_back)
      .def("__len__", &BBoxMultiCalculator::size)
      .def("__call__", &BBoxMultiCalculator::operator())
      ;

    double (PartialityCalculator::*calculate_partiality_single)(
      vec3 <double>, double, int6) const =
        &PartialityCalculator::operator();
    af::shared<double> (PartialityCalculator::*calculate_partiality_array) (
      const af::const_ref< vec3<double> >&,
      const af::const_ref<double> &,
      const af::const_ref<int6>&) const =
        &PartialityCalculator::operator();

    class_ <PartialityCalculator> ("PartialityCalculator", no_init)
      .def(init <const Beam&,
                 const Goniometer&,
                 const Scan&,
                 double > ((
        arg("beam"),
        arg("goniometer"),
        arg("scan"),
        arg("delta_mosaicity"))))
      .def("__call__", calculate_partiality_single, (
        arg("s1"), arg("phi"), arg("bbox")))
      .def("__call__", calculate_partiality_array, (
        arg("s1"), arg("phi"), arg("bbox")));

    class_ <PartialityMultiCalculator>("PartialityMultiCalculator")
      .def("append", &PartialityMultiCalculator::push_back)
      .def("__len__", &PartialityMultiCalculator::size)
      .def("__call__", &PartialityMultiCalculator::operator())
      ;
  }

}}}} // namespace = dials::algorithms::shoebox::boost_python
