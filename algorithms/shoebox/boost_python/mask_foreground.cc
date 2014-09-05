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
#include <dials/algorithms/shoebox/mask_foreground.h>

namespace dials { namespace algorithms { namespace shoebox {
  namespace boost_python {

  using namespace boost::python;

  void export_mask_foreground()
  {

    class_ <MaskForeground> ("MaskForeground", no_init)
      .def(init<const Beam&, const Detector&,
                const Goniometer&, const Scan&,
                const af::const_ref<double>&,
                const af::const_ref<double>&>((
        arg("beam"), arg("detector"),
        arg("gonio"), arg("scan"),
        arg("delta_d"), arg("delta_m"))))
      .def(init<const Beam&, const Detector&,
                const Goniometer&, const Scan&,
                double, double>((
        arg("beam"), arg("detector"),
        arg("gonio"), arg("scan"),
        arg("delta_d"), arg("delta_m"))))
      .def("__call__", &MaskForeground::single)
      .def("__call__", &MaskForeground::array);

    class_ <MaskMultiForeground>("MaskMultiForeground")
      .def("append", &MaskMultiForeground::push_back)
      .def("__len__", &MaskMultiForeground::size)
      .def("__call__", &MaskMultiForeground::operator())
      ;
  }

}}}} // namespace = dials::algorithms::shoebox::boost_python
