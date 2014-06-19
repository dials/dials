/*
 * plane.cc
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
#include <dials/algorithms/background/plane.h>

namespace dials { namespace algorithms { namespace background {
  namespace boost_python {

  using namespace boost::python;

  void export_plane()
  {
    class_<PlaneModel>("PlaneModel", no_init)
      .def(init<
          const af::const_ref< int, af::c_grid<2> > &,
          af::ref< int, af::c_grid<2> >,
          double,
          double>((
        arg("data"),
        arg("mask"),
        arg("fraction"),
        arg("nsigma"))))
      .def("a", &PlaneModel::a)
      .def("b", &PlaneModel::b)
      .def("c", &PlaneModel::c)
      .def("rmsd", &PlaneModel::rmsd);
  }

}}}} // namespace = dials::algorithms::background::boost_python

