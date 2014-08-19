/*
 * preprocessor.cc
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
#include <dials/algorithms/integration/preprocessor.h>

using namespace boost::python;

namespace dials { namespace algorithms { namespace boost_python {

  void export_preprocessor() {

    class_<PowderRingFilter>("PowderRingFilter", no_init)
      .def(init< cctbx::uctbx::unit_cell,
                 cctbx::sgtbx::space_group,
                 double,
                 double >((
        arg("unit_cell"),
        arg("space_group"),
        arg("d_min"),
        arg("width"))))
      .def("d_spacings", &PowderRingFilter::d_spacings)
      .def("d_min", &PowderRingFilter::d_min)
      .def("width", &PowderRingFilter::width)
      ;

    class_<Preprocessor>("Preprocessor", no_init)
      .def(init<af::reflection_table>())
      .def("summary", &Preprocessor::summary)
      ;
  }

}}} // namespace = dials::algorithms::boost_python

