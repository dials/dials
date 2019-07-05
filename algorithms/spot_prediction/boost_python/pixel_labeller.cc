/*
 * pixel_labeller.cc
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
#include <dials/algorithms/spot_prediction/pixel_labeller.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  static af::shared<cctbx::miller::index<> > label(const PixelLabeller &self,
                                                   mat3<double> A,
                                                   std::size_t panel_number) {
    af::c_grid<2> size = self.panel_size(panel_number);
    af::shared<cctbx::miller::index<> > result(size[0] * size[1]);
    self.label(result.ref(), A, panel_number);
    return result;
  }

  void export_pixel_labeller() {
    class_<PixelLabeller>("PixelLabeller", no_init)
      .def(init<BeamBase &, Detector>())
      .def("label", &PixelLabeller::label)
      .def("label", label);
    ;
  }

}}}  // namespace dials::algorithms::boost_python
