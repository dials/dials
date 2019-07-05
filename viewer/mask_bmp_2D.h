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
//#define DST_BTW_LIN 14
#define DST_BTW_LIN 43
#include <iostream>
#include <string>
#include <scitbx/array_family/flex_types.h>

using scitbx::af::flex_double;
using scitbx::af::flex_grid;
using scitbx::af::flex_int;

int get_mask_img_array(int (&mask_bw_img)[PX_SCALE][PX_SCALE][4]) {
  int err_cod = 0;

  // cleaning mask and painting borders
  for (int dpt = 0; dpt < 4; dpt++) {
    for (int row = 0; row < PX_SCALE; row++) {
      for (int col = 0; col < PX_SCALE; col++) {
        mask_bw_img[col][row][dpt] = 0;
      }
    }
  }

  // painting borders
  for (int dpt = 0; dpt < 4; dpt++) {
    for (int pos = 0; pos < PX_SCALE; pos++) {
      mask_bw_img[0][pos][dpt] = 1;
      mask_bw_img[PX_SCALE - 1][pos][dpt] = 1;
      mask_bw_img[pos][0][dpt] = 1;
      mask_bw_img[pos][PX_SCALE - 1][dpt] = 1;
    }
  }

  // painting diagonal lines from left top to right bottom
  for (int row = 0; row < PX_SCALE; row++) {
    for (int col = 0; col < PX_SCALE; col++) {
      for (int dg_pos = -84; dg_pos < 85; dg_pos += DST_BTW_LIN) {
        if (row == dg_pos + col) {
          mask_bw_img[col][row][0] = 1;
        }
      }
    }
  }

  // painting diagonal lines from left bottom to right top
  for (int row = 0; row < PX_SCALE; row++) {
    for (int col = 0; col < PX_SCALE; col++) {
      for (int dg_pos = 0; dg_pos < 175; dg_pos += DST_BTW_LIN) {
        if (row == dg_pos - col) {
          mask_bw_img[col][row][1] = 1;
        }
      }
    }
  }

  // painting horizontal lines
  for (int row = 0; row < PX_SCALE; row++) {
    for (int col = 0; col < PX_SCALE; col += DST_BTW_LIN) {
      mask_bw_img[col][row][2] = 1;
    }
  }

  // painting vertical lines as Background static mask
  for (int row = 0; row < PX_SCALE; row += DST_BTW_LIN) {
    for (int col = 0; col < PX_SCALE; col++) {
      mask_bw_img[col][row][3] = 1;
    }
  }

  return err_cod;
}

#endif
