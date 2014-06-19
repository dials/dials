/*
 * integrator_2d.cc
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
#include <dials/algorithms/integration/integrator_2d.h>

using namespace boost::python;

namespace dials { namespace algorithms { namespace boost_python {

  void export_integrator_2d() {

    class_<Integrator2DSpec>("Integrator2DSpec")
      ;

    class_<Integrator2D>("Integrator2D", no_init)
      .def(init<const Integrator2DSpec&>())
      ;
  }

}}} // namespace = dials::algorithms::boost_python
