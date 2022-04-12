/*
 * RGB_2d_ext.cc
 *
 *  Copyright (C) 2015 Diamond Light Source
 *
 *  Author: Luis Fuentes-Montero (Luiso)
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>

namespace dials { namespace viewer { namespace boost_python {
  using namespace boost::python;
  void export_dials_viewer();
  BOOST_PYTHON_MODULE(dials_viewer_ext) {
    export_dials_viewer();
  }
}}}  // namespace dials::viewer::boost_python
