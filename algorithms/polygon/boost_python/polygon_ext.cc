/*
 * polygon_ext.cc
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
#include <scitbx/vec2.h>
#include <dials/algorithms/polygon/area.h>
#include <dials/array_family/scitbx_shared_and_versa.h>

namespace dials { namespace algorithms { namespace polygon { namespace boost_python {

  using namespace boost::python;

  using scitbx::vec2;

  void export_area() {
    def("simple_area", &simple_area<af::shared<vec2<double> > >, (arg("poly")));
  }

  BOOST_PYTHON_MODULE(dials_algorithms_polygon_ext) {
    export_area();
  }

}}}}  // namespace dials::algorithms::polygon::boost_python
