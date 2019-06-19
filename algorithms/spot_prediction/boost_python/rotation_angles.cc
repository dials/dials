/*
 * rotation_angles.cc
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
#include <dials/algorithms/spot_prediction/rotation_angles.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_rotation_angles() {
    scitbx::vec2<double> (RotationAngles::*calculate_pstar0)(scitbx::vec3<double>)
      const = &RotationAngles::operator();

    scitbx::vec2<double> (RotationAngles::*calculate_miller)(
      cctbx::miller::index<>, scitbx::mat3<double>) const = &RotationAngles::operator();

    class_<RotationAngles>("RotationAngles", no_init)
      .def(init<scitbx::vec3<double>, scitbx::vec3<double> >(
        (arg("beam_direction"), arg("rotation_axis"))))
      .def("__call__", calculate_pstar0)
      .def("__call__", calculate_miller);
  }

}}}  // namespace dials::algorithms::boost_python
