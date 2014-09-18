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
#include <dials/algorithms/profile_model/gaussian_rs/bbox_calculator.h>

namespace dials { namespace algorithms { namespace shoebox {
  namespace boost_python {

  using namespace boost::python;

  void export_bbox_calculator()
  {
    class_ <BBoxCalculatorIface, boost::noncopyable>(
        "BBoxCalculatorIface", no_init)
      .def("__call__", &BBoxCalculatorIface::single, (
        arg("s1"),
        arg("frame"),
        arg("panel")))
      .def("__call__", &BBoxCalculatorIface::array, (
        arg("s1"),
        arg("frame"),
        arg("panel")))
      ;

    class_ <BBoxCalculator3D, bases<BBoxCalculatorIface> >(
        "BBoxCalculator3D", no_init)
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
      ;

    class_ <BBoxCalculator2D, bases<BBoxCalculatorIface> >(
        "BBoxCalculator2D", no_init)
      .def(init <const Beam&,
                 const Detector&,
                 double,
                 double > ((
        arg("beam"),
        arg("detector"),
        arg("delta_divergence"),
        arg("delta_mosaicity"))))
      ;


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
