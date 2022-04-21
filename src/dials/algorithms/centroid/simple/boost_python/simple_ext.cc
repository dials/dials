/*
 * simple_ext.cc
 *
 *  copyright (c) 2013 diamond light source
 *
 *  author: james parkhurst
 *
 *  this code is distributed under the bsd license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/algorithms/centroid/simple/algorithm.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  static void add_detector(Centroider &self, const Detector &detector) {
    self.add(detector);
  }

  static void add_detector_and_scan(Centroider &self,
                                    const Detector &detector,
                                    const Scan &scan) {
    self.add(detector, scan);
  }

  BOOST_PYTHON_MODULE(dials_algorithms_centroid_simple_ext) {
    class_<Centroider>("Centroider")
      .def("add", &add_detector)
      .def("add", &add_detector_and_scan)
      .def("__call__", &Centroider::shoebox)
      .def("__call__", &Centroider::volume<float>);
  }

}}}  // namespace dials::algorithms::boost_python
