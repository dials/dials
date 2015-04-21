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
    double max = 1, min = -1, loc_cel, dif = 0;

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

    dif = max - min;

    for (int row = 0; row < nrow ; row++) {
      for (int col = 0; col < ncol; col++) {
        bmp_dat(row, col) = 255.0 * ( ( bmp_dat(row, col) - min ) / dif );
      }
    }


    /*

  width = np.size( data2d_scale[0:1, :] )
  height = np.size( data2d_scale[:, 0:1] )
  img_array = np.empty( (height ,width, 3), 'int')

  img_array_r = np.empty( (height, width), 'int')
  img_array_g = np.empty( (height, width), 'int')
  img_array_b = np.empty( (height, width), 'int')

  scaled_i = np.empty( (height, width), 'int')

  red_byte = np.empty( (255 * 3), 'int')
  green_byte = np.empty( (255 * 3), 'int')
  blue_byte = np.empty( (255 * 3), 'int')

  for i in xrange(255):
    red_byte[i] = i
    green_byte[i + 255] = i
    blue_byte[i + 255 * 2 ] = i


  red_byte[255:255 * 3] = 255
  green_byte[0:255] = 0
  green_byte[255 * 2 : 255 * 3] = 255
  blue_byte[0:255 * 2] = 0

  blue_byte[764] = 255
  red_byte[764] = 255
  green_byte[764] = 255

  scaled_i[:,:] = data2d_scale[:,:]

  img_array_r[:,:] = scaled_i[:,:]
  for x in np.nditer(img_array_r[:,:], op_flags=['readwrite'], flags=['external_loop']):
    x[...] = red_byte[x]

  img_array_g[:,:] = scaled_i[:,:]
  for x in np.nditer(img_array_g[:,:], op_flags=['readwrite'], flags=['external_loop']):
    x[...] = green_byte[x]

  img_array_b[:,:] = scaled_i[:,:]
  for x in np.nditer(img_array_b[:,:], op_flags=['readwrite'], flags=['external_loop']):
    x[...] = blue_byte[x]

  img_array[:, :, 0] = img_array_r[:,:]
  img_array[:, :, 1] = img_array_g[:,:]
  img_array[:, :, 2] = img_array_b[:,:]
  print ("From Python call")
    */

  return bmp_dat;
  }

}}}

#endif
