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

  void export_pixel_to_miller_index()
  {
    class_ <PixelToMillerIndex> ("PixelToMillerIndex", no_init)
      .def(init <
          const Beam&,
          const Detector&,
          const Goniometer&,
          const Scan&,
          const Crystal&>())
      .def("h", &PixelToMillerIndex::h)
      ;
  }

}}} // namespace = dials::algorithms::boost_python
