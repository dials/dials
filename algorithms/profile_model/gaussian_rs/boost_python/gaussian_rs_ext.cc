/*
 * gaussian_rs_ext.cc
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

  void export_bbox_calculator();
  void export_mask_foreground();
  void export_mask_foreground_2d();

  BOOST_PYTHON_MODULE(dials_algorithms_profile_model_gaussian_rs_ext)
  {
    export_bbox_calculator();
    export_mask_foreground();
    export_mask_foreground_2d();
  }

}}}} // namespace = dials::algorithms::shoebox::boost_python
