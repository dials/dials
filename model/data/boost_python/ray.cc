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

  void export_ray()
  {
    class_<Ray>("Ray")
      .def_readwrite("s1", &Ray::s1)
      .def_readwrite("angle", &Ray::angle)
      .def_readwrite("entering", &Ray::entering);
  }

}}} // namespace dials::model::boost_python
