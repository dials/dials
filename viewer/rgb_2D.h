/*
 * rgb_2d.h
 *
 *  Copyright (C) 2015 Diamond Light Source
 *
 *  Author: Luis Fuentes-Montero (Luiso)
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef DIALS_RGB_IMG_BUILDER_H
#define DIALS_RGB_IMG_BUILDER_H
#include <iostream>
#include <scitbx/array_family/flex_types.h>

namespace dials { namespace viewer { namespace boost_python {
  using scitbx::af::flex_double;
  using scitbx::af::flex_grid;


  flex_double gen_img(flex_double & data2d) {
    int ncol=data2d.accessor().all()[1];
    int nrow=data2d.accessor().all()[0];
    double max = 1, min = -1, loc_cel;

    flex_double bmp_dat(flex_grid<>(nrow, ncol),0);

    for (int row = 0; row < nrow ; row++) {
      for (int col = 0; col < ncol; col++) {

        loc_cel = data2d(row, col);
        if(row == 0 and col == 0){
          max = loc_cel;
          min = loc_cel;
        } else {
          if(loc_cel > max){
            max = loc_cel;
            }
          if(loc_cel < min){
            min = loc_cel;
            }
          }
        bmp_dat(row, col) = data2d(row, col);
      }


    }

    std::cout << "\n max = "<< max << "\n";
    std::cout << "\n min = "<< min << "\n \n";

  return bmp_dat;
  }

}}}

#endif
