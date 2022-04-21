/*
 * corrections.cc
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
#include <dials/algorithms/integration/corrections.h>

using namespace boost::python;

namespace dials { namespace algorithms { namespace boost_python {

  void export_corrections() {
    def("lp_correction",
        &lp_correction,
        (arg("s0"), arg("pn"), arg("pf"), arg("m2"), arg("s1")));

    def("lp_correction",
        &stills_lp_correction,
        (arg("s0"), arg("pn"), arg("pf"), arg("s1")));

    def("qe_correction", &qe_correction, (arg("mu"), arg("t0"), arg("s1"), arg("n")));

    class_<Corrections>("Corrections", no_init)
      .def(init<const BeamBase&, const Goniometer&, const Detector&>())
      .def(init<const BeamBase&, const Detector&>())
      .def("lp", &Corrections::lp, (arg("s1")))
      .def("qe", &Corrections::qe, (arg("s1"), arg("panel")));

    class_<CorrectionsMulti>("CorrectionsMulti")
      .def("append", &CorrectionsMulti::push_back)
      .def("__len__", &CorrectionsMulti::size)
      .def("lp", &CorrectionsMulti::lp)
      .def("qe", &CorrectionsMulti::qe);
  }

}}}  // namespace dials::algorithms::boost_python
