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
#include <dials/algorithms/profile_model/gaussian_rs/partiality_calculator.h>

namespace dials { namespace algorithms { namespace shoebox {
  namespace boost_python {

  using namespace boost::python;

  void export_partiality_calculator()
  {
    class_ <PartialityCalculatorIface, boost::noncopyable>(
        "PartialityCalculatorIface", no_init)
      .def("__call__", &PartialityCalculatorIface::single, (
        arg("s1"),
        arg("frame"),
        arg("bbox")))
      .def("__call__", &PartialityCalculatorIface::array, (
        arg("s1"),
        arg("frame"),
        arg("bbox")))
      ;

    class_ <PartialityCalculator3D, bases<PartialityCalculatorIface> >(
        "PartialityCalculator3D", no_init)
      .def(init <const Beam&,
                 const Goniometer&,
                 const Scan&,
                 double > ((
        arg("beam"),
        arg("goniometer"),
        arg("scan"),
        arg("delta_mosaicity"))))
      ;

    class_ <PartialityCalculator2D, bases<PartialityCalculatorIface> >(
        "PartialityCalculator2D", no_init)
      .def(init <const Beam&,
                 double > ((
        arg("beam"),
        arg("delta_mosaicity"))))
      ;


    class_ <PartialityMultiCalculator>("PartialityMultiCalculator")
      .def("append", &PartialityMultiCalculator::push_back)
      .def("__len__", &PartialityMultiCalculator::size)
      .def("__call__", &PartialityMultiCalculator::operator())
      ;
  }

}}}} // namespace = dials::algorithms::shoebox::boost_python
