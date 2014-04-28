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

    def("sigma_2d", &sigma_2d, (arg("intensity")),
        (arg("mask2d")), (arg("background2d")));

    def("add_2d", &add_2d, arg("descriptor"), arg("data2d"), arg("tmp_total"));

    def("mask_2d_interpolate", &mask_2d_interpolate,
            arg("descriptor"), arg("mask2d_in"), arg("mask2d_tmp_total"));

    def("simple_2d_add", &simple_2d_add, arg("in_data2d_one"),
        arg("in_data2d_two"));

    def("mask_add_2d", &mask_add_2d, arg("mask2d_one"),arg("mask2d_one"));


    def("fitting_2d_multile_var_build_mat", &fitting_2d_multile_var_build_mat,
        (arg("descriptor"), arg("data2d"), arg("background2d"),
         arg("profile2d"), arg("tmp_scale") ,  arg("mat_a"), arg("vec_b")));

    def("fitting_2d_partials", &fitting_2d_partials,
        (arg("descriptor"), arg("data2dmov"), arg("backg2dmov"), arg("profile2d"),
         arg("interpolation_mask2d"),      arg("sum_its") ));

    def("subtrac_bkg_2d", &subtrac_bkg_2d, (arg("data2d"), arg("background2d")));


  }

}}}
