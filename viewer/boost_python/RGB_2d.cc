/*
 * RGB_2d.cc
 *
 *  Copyright (C) 2015 Diamond Light Source
 *
 *  Author: Luis Fuentes-Montero (Luiso)
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/viewer/rgb_2D.h>

namespace dials { namespace viewer { namespace boost_python {
  using namespace boost::python;
  void export_dials_viewer() {
    // def("gen_font_img", &gen_font_img, (arg("data2d")));
    // def("gen_str_tst", &gen_str_tst, (arg("data_num")));

    class_<rgb_img>("rgb_img")

      .def("set_min_max", &rgb_img::set_min_max, (arg("new_min"), arg("new_max")))

      .def("gen_bmp",
           &rgb_img::gen_bmp,
           (arg("data2d"), arg("mask2d"), arg("show_nums"), arg("palette_num")));

    // def("tst_ref_prod", &tst_ref_prod, arg("matr01"), arg("matr02"));
  }
}}}  // namespace dials::viewer::boost_python
