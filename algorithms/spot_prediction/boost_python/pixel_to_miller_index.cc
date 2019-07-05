/*
 * pixel_to_miller_index.cc
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
#include <dials/algorithms/spot_prediction/pixel_to_miller_index.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  static vec3<double> PixelToMillerIndex_h_rotation(const PixelToMillerIndex &self,
                                                    std::size_t panel,
                                                    double x,
                                                    double y,
                                                    double z) {
    return self.h(panel, x, y, z);
  }

  static vec3<double> PixelToMillerIndex_h_stills(const PixelToMillerIndex &self,
                                                  std::size_t panel,
                                                  double x,
                                                  double y) {
    return self.h(panel, x, y);
  }

  void export_pixel_to_miller_index() {
    class_<PixelToMillerIndex>("PixelToMillerIndex", no_init)
      .def(init<const BeamBase &,
                const Detector &,
                const Goniometer &,
                const Scan &,
                const CrystalBase &>())
      .def(init<const BeamBase &, const Detector &, const CrystalBase &>())
      .def("h", &PixelToMillerIndex_h_rotation)
      .def("h", &PixelToMillerIndex_h_stills)
      .def("q", &PixelToMillerIndex::q);
  }

}}}  // namespace dials::algorithms::boost_python
