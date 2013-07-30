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

namespace dials { namespace algorithms { namespace shoebox {
  namespace boost_python {

  using namespace boost::python;

  void export_mask_code();
  void export_bbox_calculator();
  void export_find_overlapping();
  void export_mask_foreground();
  void export_masker();
  void export_populator();

  BOOST_PYTHON_MODULE(dials_algorithms_shoebox_ext)
  {
    export_mask_code();
    export_bbox_calculator();
    export_find_overlapping();
    export_mask_foreground();
    export_masker();
    export_populator();
  }

}}}} // namespace = dials::algorithms::shoebox::boost_python
