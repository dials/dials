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
#include <dials/algorithms/shoebox/mask_foreground_2d.h>

namespace dials { namespace algorithms { namespace shoebox {
  namespace boost_python {

  using namespace boost::python;

  void export_mask_foreground_2d()
  {
    class_ <MaskForeground2d> ("MaskForeground2d", no_init)
      .def(init<const Beam&,
                const Detector&,
                double>((
        arg("beam"),
        arg("detector"),
        arg("delta_d"))))
      .def("__call__", &MaskForeground2d::operator());
  }

}}}} // namespace = dials::algorithms::shoebox::boost_python
