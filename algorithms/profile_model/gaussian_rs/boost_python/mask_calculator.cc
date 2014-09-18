/*
 * mask_foreground.cc
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
#include <dials/algorithms/profile_model/gaussian_rs/mask_calculator.h>

namespace dials { namespace algorithms { namespace shoebox {
  namespace boost_python {

  using namespace boost::python;

  void export_mask_calculator()
  {
    class_ <MaskCalculatorIface, boost::noncopyable>(
        "MaskCalculatorIface", no_init)
      .def("__call__", &MaskCalculatorIface::single, (
        arg("shoebox"),
        arg("s1"),
        arg("frame"),
        arg("panel")))
      .def("__call__", &MaskCalculatorIface::array, (
        arg("shoebox"),
        arg("s1"),
        arg("frame"),
        arg("panel")))
      ;

    class_ <MaskCalculator3D, bases<MaskCalculatorIface> >(
        "MaskCalculator3D", no_init)
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

    class_ <MaskCalculator2D, bases<MaskCalculatorIface> >(
        "MaskCalculator2D", no_init)
      .def(init <const Beam&,
                 const Detector&,
                 double,
                 double > ((
        arg("beam"),
        arg("detector"),
        arg("delta_divergence"),
        arg("delta_mosaicity"))))
      ;

    class_ <MaskMultiCalculator>("MaskMultiCalculator")
      .def("append", &MaskMultiCalculator::push_back)
      .def("__len__", &MaskMultiCalculator::size)
      .def("__call__", &MaskMultiCalculator::operator())
      ;
  }

}}}} // namespace = dials::algorithms::shoebox::boost_python
