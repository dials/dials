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
  using scitbx::af::flex_grid;

  MaskBadPixels* make_mask_bad_pixels(af::versa< bool, flex_grid<> > mask) {
    return new MaskBadPixels(af::versa< bool, af::c_grid<2> >(
      mask.handle(), af::c_grid<2>(mask.accessor())));
  }

  void export_mask_bad_pixels()
  {
    class_ <MaskBadPixels> ("MaskBadPixels", no_init)
      .def("__init__", make_constructor(&make_mask_bad_pixels))
      .def("__call__", 
        (void(MaskBadPixels::*)(Reflection&)const)
          &MaskBadPixels::operator(), (
          arg("reflection")))
      .def("__call__", 
        (void(MaskBadPixels::*)(af::ref<Reflection>)const)
          &MaskBadPixels::operator(), (
          arg("reflection_list")));
            
  }

}}}} // namespace = dials::algorithms::shoebox::boost_python
