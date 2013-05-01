/*
 * peak_finding_ext.cc
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

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_label_pixels();
  void export_lui_find_peak_helper();
  void export_lui_integrate_helper();
  void export_lui_find_peak_smoothing();

  BOOST_PYTHON_MODULE(dials_algorithms_peak_finding_ext)
  {
    export_label_pixels();
    export_lui_find_peak_helper();
    export_lui_integrate_helper();
    export_lui_find_peak_smoothing();
  }

}}} // namespace = dials::algorithms::boost_python
