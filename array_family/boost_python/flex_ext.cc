/*
 * flex_ext.cc
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
#include <string>
#include <scitbx/boost_python/container_conversions.h>
#include <scitbx/array_family/boost_python/c_grid_flex_conversions.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/model/data/ray.h>
#include <dials/config.h>

namespace dials { namespace af { namespace boost_python {

  using namespace boost::python;
  using namespace scitbx::boost_python::container_conversions;

  void export_flex_int6();
  void export_flex_shoebox();
  void export_flex_centroid();
  void export_flex_intensity();
  void export_flex_observation();
  void export_flex_prediction();
  void export_flex_reflection_table();
  void export_flex_unit_cell();
  void export_flex_shoebox_extractor();
  void export_flex_binner();

  template <typename FloatType>
  std::string get_real_type();

  template <>
  std::string get_real_type<float>() {
    return "float";
  }

  template <>
  std::string get_real_type<double>() {
    return "double";
  }

  BOOST_PYTHON_MODULE(dials_array_family_flex_ext) {
    export_flex_int6();
    export_flex_shoebox();
    export_flex_centroid();
    export_flex_intensity();
    export_flex_observation();
    export_flex_prediction();
    export_flex_reflection_table();
    export_flex_unit_cell();
    export_flex_shoebox_extractor();
    export_flex_binner();

    def("get_real_type", &get_real_type<ProfileFloatType>);

    scitbx::af::boost_python::c_grid_flex_conversions<bool, af::c_grid<4> >();
    scitbx::af::boost_python::c_grid_flex_conversions<double, af::c_grid<4> >();

    tuple_mapping_fixed_capacity<af::small<model::Ray, 2> >();
  }

}}}  // namespace dials::af::boost_python
