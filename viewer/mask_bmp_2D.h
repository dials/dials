/*
 * mask_bmp_2D.h
 *
 *  Copyright (C) 2015 Diamond Light Source
 *
 *  Author: Luis Fuentes-Montero (Luiso)
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef DIALS_MASK_LOW_LEVEL_H
#define DIALS_MASK_LOW_LEVEL_H
#define PX_SCALE 85
#include <iostream>
#include <string>
#include <cmath>
#include <scitbx/array_family/flex_types.h>

using scitbx::af::flex_double;
using scitbx::af::flex_int;
using scitbx::af::flex_grid;

int get_mask_img_array( int (&mask_bw_img)[PX_SCALE][PX_SCALE][5]){
  int err_cod = 0;







  std::cout << "\n Hi from mask_bmp_2D\n";

  return err_cod;
}


#endif
