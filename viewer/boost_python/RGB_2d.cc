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
  void dials_viewer(){
    def("gen_img", &gen_img, (arg("data2d")));
  }
}}}
