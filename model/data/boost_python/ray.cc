/*
 * ray.cc
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
#include <scitbx/vec3.h>
#include <dials/model/data/ray.h>

namespace dials { namespace model { namespace boost_python {

  using namespace boost::python;
  using scitbx::vec2;
  using scitbx::vec3;

  static vec3<double> get_s1(const Ray &ray) {
    return ray.s1;
  }

  static void set_s1(Ray &ray, const vec3<double> &s1) {
    ray.s1 = s1;
  }

  void export_ray() {
    class_<Ray>("Ray")
      .add_property("s1", get_s1, set_s1)
      .def_readwrite("angle", &Ray::angle)
      .def_readwrite("entering", &Ray::entering);
  }

}}}  // namespace dials::model::boost_python
