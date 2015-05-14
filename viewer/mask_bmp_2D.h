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

  for(int dpt = 0; dpt < 5; dpt++){
    for(int row = 0; row < PX_SCALE; row++){
      for(int col = 0; col < PX_SCALE; col++){
        mask_bw_img[col][row][dpt] = 0;
      }
    }
  }

  for(int row = 0; row < PX_SCALE; row++){
    for(int col = 0; col < PX_SCALE; col++){
      for(int dg_pos = -80; dg_pos < 85; dg_pos+= 8){
        if(row == col + dg_pos){
          mask_bw_img[col][row][0] = 1;
        }
      }

    }
  }

  for(int row = 0; row < PX_SCALE; row++){
    for(int col = 0; col < PX_SCALE; col++){
      if(row == 2 * col){
        mask_bw_img[col][row][1] = 1;
      }else{
        mask_bw_img[col][row][1] = 0;
      }
    }
  }

  for(int row = 0; row < PX_SCALE; row++){
    for(int col = 0; col < PX_SCALE; col++){
      if(2 * row == col){
        mask_bw_img[col][row][2] = 1;
      }else{
        mask_bw_img[col][row][2] = 0;
      }
    }
  }

  for(int row = 0; row < PX_SCALE; row++){
    for(int col = 0; col < PX_SCALE; col++){
      if(row == col + 5){
        mask_bw_img[col][row][3] = 1;
      }else{
        mask_bw_img[col][row][3] = 0;
      }
    }
  }
  for(int row = 0; row < PX_SCALE; row++){
    for(int col = 0; col < PX_SCALE; col++){
      if(row == col - 5){
        mask_bw_img[col][row][4] = 1;
      }else{
        mask_bw_img[col][row][4] = 0;
      }
    }
  }

  std::cout << "\n Hi from mask_bmp_2D\n";

  return err_cod;
}


#endif
