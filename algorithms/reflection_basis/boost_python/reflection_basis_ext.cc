/*
 * reflection_basis_ext.cc
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
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/algorithms/reflection_basis/coordinate_system.h>

namespace dials { namespace algorithms { namespace reflection_basis {
  namespace boost_python {

  using namespace boost::python;

  /**
   * Helper function to calculate zeta for an array of s1.
   */
  af::shared<double> zeta_factor_array(
      vec3<double> m2, vec3<double> s0,
      const af::const_ref< vec3<double> > &s1) {
    af::shared<double> result(s1.size());
    for (std::size_t i = 0; i < s1.size(); ++i) {
      result[i] = zeta_factor(m2, s0, s1[i]);
    }
    return result;
  }

  void export_coordinate_system()
  {
    // Export zeta factor functions
    def("zeta_factor",
      (double (*)(vec3<double>, vec3<double>, vec3<double>))&zeta_factor,
      (arg("m2"), arg("s0"), arg("s1")));
    def("zeta_factor",
      (double (*)(vec3<double>, vec3<double>))&zeta_factor,
      (arg("m2"), arg("e1")));
    def("zeta_factor", &zeta_factor_array, (
      arg("m2"), arg("s0"), arg("s1")));

    // Export coordinate system
    class_<CoordinateSystem>("CoordinateSystem", no_init)
      .def(init<vec3<double>,
                vec3<double>,
                vec3<double>,
                double>((
        arg("m2"),
        arg("s0"),
        arg("s1"),
        arg("phi"))))
      .def("m2", &CoordinateSystem::m2)
      .def("s0", &CoordinateSystem::s0)
      .def("s1", &CoordinateSystem::s1)
      .def("phi", &CoordinateSystem::phi)
      .def("p_star", &CoordinateSystem::p_star)
      .def("e1_axis", &CoordinateSystem::e1_axis)
      .def("e2_axis", &CoordinateSystem::e2_axis)
      .def("e3_axis", &CoordinateSystem::e3_axis)
      .def("zeta", &CoordinateSystem::zeta)
      .def("lorentz_inv", &CoordinateSystem::lorentz_inv)
      .def("lorentz", &CoordinateSystem::lorentz)
      .def("path_length_increase", &CoordinateSystem::path_length_increase)
      .def("limits", &CoordinateSystem::limits)
      .def("from_beam_vector",
        &CoordinateSystem::from_beam_vector)
      .def("from_rotation_angle",
        &CoordinateSystem::from_rotation_angle)
      .def("from_rotation_angle_fast",
        &CoordinateSystem::from_rotation_angle_fast)
      .def("from_beam_vector_and_rotation_angle",
        &CoordinateSystem::from_beam_vector_and_rotation_angle)
      .def("to_beam_vector",
        &CoordinateSystem::to_beam_vector)
      .def("to_rotation_angle",
        &CoordinateSystem::to_rotation_angle)
      .def("to_rotation_angle_fast",
        &CoordinateSystem::to_rotation_angle_fast)
      .def("to_beam_vector_and_rotation_angle",
        &CoordinateSystem::to_beam_vector_and_rotation_angle);
  }

  BOOST_PYTHON_MODULE(dials_algorithms_reflection_basis_ext)
  {
    export_coordinate_system();
  }

}}}} // namespace = dials::algorithms::reflection_basis::boost_python
