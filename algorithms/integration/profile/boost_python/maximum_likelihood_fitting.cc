/*
 * maximum_likelihood_fitting.cc
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
#include <boost/python/iterator.hpp>
#include <dials/algorithms/integration/profile/maximum_likelihood_fitting.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_maximum_likelihood_fitting() {
    class_<MLPoisson2>("MLPoisson2", no_init)
      .def(init< const af::const_ref<double>&,
                 const af::const_ref<double>&,
                 const af::const_ref<double>&,
                 double,
                 double,
                 double >())
      .def("step", &MLPoisson2::step)
      .def("solve", &MLPoisson2::solve)
      .def("S1", &MLPoisson2::S1)
      .def("S2", &MLPoisson2::S2)
      ;

    class_<MLPoisson2Stepper>("MLPoisson2Stepper", no_init)
      .def(init< const af::const_ref<double>&,
                 const af::const_ref<double>&,
                 const af::const_ref<double>&,
                 vec2<double> >())
      .def("step", &MLPoisson2Stepper::step)
      .def("X", &MLPoisson2Stepper::X)
      .def("B", &MLPoisson2Stepper::B);
  }

}}} // namespace = dials::algorithms::boost_python

