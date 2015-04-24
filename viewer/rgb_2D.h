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
  using scitbx::af::flex_int;
  using scitbx::af::flex_grid;

  /*
  flex_int gen_img(flex_double & data2d) {

    int ncol=data2d.accessor().all()[1];
    int nrow=data2d.accessor().all()[0];
    flex_int bmp_dat(flex_grid<>(nrow, ncol, 3),0);

    return bmp_dat;
  }
  */

  class rgb_img
  {

    private:

      int red_byte[255 * 3];
      int green_byte[255 * 3];
      int blue_byte[255 * 3];
    public:
      rgb_img() {

        for (int i = 0; i < 255; i++){
          red_byte[i] = i;
          green_byte[i + 255] = i;
          blue_byte[i + 255 * 2 ] = i;
        }

        for (int i = 255; i < 255 * 3; i++){
          red_byte[i] = 255;
        }

        for (int i = 0; i < 255; i++){
          green_byte[i] = 0;
        }

        for (int i = 255 * 2; i < 255 * 3; i++){
          green_byte[i] = 255;
        }

        for (int i = 0; i < 255 * 2; i++){
          blue_byte[i] = 0;
        }

        blue_byte[764] = 255;
        green_byte[764] = 255;
        red_byte[764] = 255;

      }

      flex_int gen_bmp(flex_double & data2d) {

        int ncol=data2d.accessor().all()[1];
        int nrow=data2d.accessor().all()[0];
        double max = 1, min = -1, loc_cel, dif = 0;

        flex_double scaled_array(flex_grid<>(nrow, ncol),0);

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
            scaled_array(row, col) = data2d(row, col);
          }
        }

        std::cout << "\n max = "<< max << "\n";
        std::cout << "\n min = "<< min << "\n \n";

        dif = max - min;

        for (int row = 0; row < nrow ; row++) {
          for (int col = 0; col < ncol; col++) {
            scaled_array(row, col) = 255.0 * 3 * ( ( scaled_array(row, col) - min ) / dif );
          }
        }


        flex_int bmp_dat(flex_grid<>(nrow, ncol, 3),0);

        for (int row = 0; row < nrow ; row++) {
          for (int col = 0; col < ncol; col++) {
            bmp_dat(row, col, 0) = red_byte[int(scaled_array(row, col))];
            bmp_dat(row, col, 1) = green_byte[int(scaled_array(row, col))];
            bmp_dat(row, col, 2) = blue_byte[int(scaled_array(row, col))];
          }
        }

      return bmp_dat;
      }

  };



}}}

#endif
