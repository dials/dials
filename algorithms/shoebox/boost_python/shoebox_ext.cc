/*
 * shoebox_ext.cc
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

  void export_bbox_calculator();
  void export_find_overlapping_reflections();
  void export_shoebox_masker();

  BOOST_PYTHON_MODULE(dials_algorithms_shoebox_ext)
  {
    export_bbox_calculator();
    export_find_overlapping_reflections();
    export_shoebox_masker();
  }

}}} // namespace = dials::algorithms::boost_python
