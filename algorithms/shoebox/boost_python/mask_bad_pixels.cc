/*
 * mask_bad_pixels.cc
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
#include <dials/algorithms/shoebox/mask_bad_pixels.h>

namespace dials { namespace algorithms { namespace shoebox {
  namespace boost_python {

  using namespace boost::python;

  void export_mask_bad_pixels()
  {
    class_ <MaskBadPixels> ("MaskBadPixels", no_init)
      .def(init<const flex_bool&>((
        arg("detector_mask"))))
      .def("__call__", 
        (void(MaskBadPixels::*)(Reflection&)const)
          &MaskBadPixels::operator(), (
          arg("reflection")))
      .def("__call__", 
        (void(MaskBadPixels::*)(ReflectionList&)const)
          &MaskBadPixels::operator(), (
          arg("reflection_list")));
            
  }

}}}} // namespace = dials::algorithms::shoebox::boost_python
