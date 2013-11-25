/*
 * lui_2d_integration.cc
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: Luis Fuentes-Montero
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/algorithms/integration/lui_2d_integration.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_luiso_s_2d_integration() {
    def("raw_2d_cut", &raw_2d_cut, (arg("data2d")), (arg("mask2d")),
            (arg("background2d")));
    def("add_2d", &add_2d, arg("descriptor"), arg("data2d"), arg("tmp_total"));

    def("fitting_2d", &fitting_2d, (arg("descriptor"), arg("data2d"), arg("background2d"), arg("profile2d") ));
    def("subtrac_bkg_2d", &subtrac_bkg_2d, (arg("data2d"), arg("background2d")));
  }

}}}
