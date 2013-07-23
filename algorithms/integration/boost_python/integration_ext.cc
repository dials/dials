/*
 * integration_ext.cc
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

  void export_helpers();
  void export_bbox_calculator();
  void export_find_overlapping_reflections();
  void export_reflection_profile_helpers();
  void export_shoebox_masker();
  void export_luiso_s_2d_integration();
  void export_summation();
  void export_summation_reciprocal_space();
  void export_profile_fitting_reciprocal_space();

  BOOST_PYTHON_MODULE(dials_algorithms_integration_ext)
  {
    export_helpers();
    export_bbox_calculator();
    export_find_overlapping_reflections();
    export_reflection_profile_helpers();
    export_shoebox_masker();
    export_luiso_s_2d_integration();
    export_summation();
    export_summation_reciprocal_space();
    export_profile_fitting_reciprocal_space();
  }

}}} // namespace = dials::algorithms::boost_python
