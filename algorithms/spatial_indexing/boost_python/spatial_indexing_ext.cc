/*
 * spatial_indexing_ext.cc
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

  void export_quadtree();
  void export_octree();

  BOOST_PYTHON_MODULE(dials_algorithms_spatial_indexing_ext) {
    export_quadtree();
    export_octree();
  }

}}}  // namespace dials::algorithms::boost_python
