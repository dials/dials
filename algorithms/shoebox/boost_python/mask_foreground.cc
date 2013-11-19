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
    void (MaskForeground::*call_w_reflection)(Reflection&)const =
      &MaskForeground::operator();

    void (MaskForeground::*call_w_reflection_list)(af::ref<Reflection>)const =
      &MaskForeground::operator();

    class_ <MaskForeground> ("MaskForeground", no_init)
      .def(init<const Beam&, const Detector&,
                const Goniometer&, const Scan&,
                double, double>((
        arg("beam"), arg("detector"),
        arg("gonio"), arg("scan"),
        arg("delta_d"), arg("delta_m"))))
      .def("__call__", call_w_reflection)
      .def("__call__", call_w_reflection_list);
  }

}}}} // namespace = dials::algorithms::shoebox::boost_python
