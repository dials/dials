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

#ifndef DIALS_RGB_LOW_LEVEL_H
#define DIALS_RGB_LOW_LEVEL_H
#include <iostream>
#include <string>
#include <cmath>
#include <scitbx/array_family/flex_types.h>


  using scitbx::af::flex_double;
  using scitbx::af::flex_int;
  using scitbx::af::flex_grid;

  int get_font_img_array( int (&font_bw_img)[7][7][16]){

    int err_cod = 0;

    int arr_2d_0[7][7] = {{0,0,0,0,0,0,0},
                          {0,0,1,1,1,1,0},
                          {0,1,1,0,0,1,1},
                          {0,1,0,0,1,0,1},
                          {0,1,0,1,0,0,1},
                          {0,0,1,1,1,1,0},
                          {0,0,0,0,0,0,0}};
    for (int row = 0; row < 7 ; row++) {
      for (int col = 0; col < 7; col++) {
        font_bw_img[col][row][0] = arr_2d_0[col][row];
      }
    }

    int arr_2d_1[7][7] = {{0,0,0,0,0,0,0},
                          {0,1,1,1,0,0,0},
                          {0,0,0,1,0,0,0},
                          {0,0,0,1,0,0,0},
                          {0,0,0,1,0,0,0},
                          {0,1,1,1,1,1,0},
                          {0,0,0,0,0,0,0}};
    for (int row = 0; row < 7 ; row++) {
      for (int col = 0; col < 7; col++) {
        font_bw_img[col][row][1] = arr_2d_1[col][row];
      }
    }

    int arr_2d_2[7][7] = {{0,0,0,0,0,0,0},
                          {0,0,1,1,1,1,0},
                          {0,1,0,0,0,1,1},
                          {0,0,0,0,1,1,0},
                          {0,0,0,1,1,0,0},
                          {0,1,1,1,1,1,1},
                          {0,0,0,0,0,0,0}};
    for (int row = 0; row < 7 ; row++) {
      for (int col = 0; col < 7; col++) {
        font_bw_img[col][row][2] = arr_2d_2[col][row];
      }
    }

    int arr_2d_3[7][7] = {{0,0,0,0,0,0,0},
                          {0,0,1,1,1,1,0},
                          {0,1,0,0,0,1,1},
                          {0,0,0,1,1,1,0},
                          {0,0,0,0,0,1,1},
                          {0,1,1,1,1,1,0},
                          {0,0,0,0,0,0,0}};
    for (int row = 0; row < 7 ; row++) {
      for (int col = 0; col < 7; col++) {
        font_bw_img[col][row][3] = arr_2d_3[col][row];
      }
    }

    int arr_2d_4[7][7] = {{0,0,0,0,0,0,0},
                          {0,0,0,0,1,1,0},
                          {0,0,0,1,1,1,0},
                          {0,0,1,1,0,1,0},
                          {0,1,1,1,1,1,1},
                          {0,0,0,0,0,1,0},
                          {0,0,0,0,0,0,0}};
    for (int row = 0; row < 7 ; row++) {
      for (int col = 0; col < 7; col++) {
        font_bw_img[col][row][4] = arr_2d_4[col][row];
      }
    }

    int arr_2d_5[7][7] = {{0,0,0,0,0,0,0},
                          {0,1,1,1,1,1,0},
                          {0,1,0,0,0,0,0},
                          {0,1,1,1,1,1,1},
                          {0,0,0,0,0,0,1},
                          {0,1,1,1,1,1,1},
                          {0,0,0,0,0,0,0}};
    for (int row = 0; row < 7 ; row++) {
      for (int col = 0; col < 7; col++) {
        font_bw_img[col][row][5] = arr_2d_5[col][row];
      }
    }

    int arr_2d_6[7][7] = {{0,0,0,0,0,0,0},
                          {0,0,1,1,1,1,0},
                          {0,1,1,0,0,0,0},
                          {0,1,1,1,1,1,1},
                          {0,1,1,0,0,0,1},
                          {0,0,1,1,1,1,1},
                          {0,0,0,0,0,0,0}};
    for (int row = 0; row < 7 ; row++) {
      for (int col = 0; col < 7; col++) {
        font_bw_img[col][row][6] = arr_2d_6[col][row];
      }
    }

    int arr_2d_7[7][7] = {{0,0,0,0,0,0,0},
                          {0,1,1,1,1,1,1},
                          {0,0,0,0,0,1,1},
                          {0,0,0,0,1,1,0},
                          {0,0,0,1,1,0,0},
                          {0,0,1,1,0,0,0},
                          {0,0,0,0,0,0,0}};
    for (int row = 0; row < 7 ; row++) {
      for (int col = 0; col < 7; col++) {
        font_bw_img[col][row][7] = arr_2d_7[col][row];
      }
    }

    int arr_2d_8[7][7] = {{0,0,0,0,0,0,0},
                          {0,0,1,1,1,1,0},
                          {0,1,1,0,0,1,1},
                          {0,0,1,1,1,1,0},
                          {0,1,1,0,0,1,1},
                          {0,0,1,1,1,1,0},
                          {0,0,0,0,0,0,0}};
    for (int row = 0; row < 7 ; row++) {
      for (int col = 0; col < 7; col++) {
        font_bw_img[col][row][8] = arr_2d_8[col][row];
      }
    }

    int arr_2d_9[7][7] = {{0,0,0,0,0,0,0},
                          {0,0,1,1,1,1,0},
                          {0,1,1,0,0,1,1},
                          {0,1,1,1,1,1,1},
                          {0,0,0,0,0,1,1},
                          {0,1,1,1,1,1,0},
                          {0,0,0,0,0,0,0}};
    for (int row = 0; row < 7 ; row++) {
      for (int col = 0; col < 7; col++) {
        font_bw_img[col][row][9] = arr_2d_9[col][row];
      }
    }

    int arr_2d_point[7][7] = {{0,0,0,0,0,0,0},
                              {0,0,0,0,0,0,0},
                              {0,0,0,0,0,0,0},
                              {0,0,0,0,0,0,0},
                              {0,0,0,1,1,0,0},
                              {0,0,0,1,1,0,0},
                              {0,0,0,0,0,0,0}};
    for (int row = 0; row < 7 ; row++) {
      for (int col = 0; col < 7; col++) {
        font_bw_img[col][row][10] = arr_2d_point[col][row];
      }
    }

    int arr_2d_e[7][7] = {{0,0,0,0,0,0,0},
                          {0,0,1,1,1,1,0},
                          {0,1,1,0,0,1,0},
                          {0,1,1,1,1,1,0},
                          {0,1,1,0,0,0,0},
                          {0,0,1,1,1,1,0},
                          {0,0,0,0,0,0,0}};
    for (int row = 0; row < 7 ; row++) {
      for (int col = 0; col < 7; col++) {
        font_bw_img[col][row][11] = arr_2d_e[col][row];
      }
    }

    int arr_2d_plus[7][7] =  {{0,0,0,0,0,0,0},
                              {0,0,0,1,0,0,0},
                              {0,0,0,1,0,0,0},
                              {0,1,1,1,1,1,0},
                              {0,0,0,1,0,0,0},
                              {0,0,0,1,0,0,0},
                              {0,0,0,0,0,0,0}};
    for (int row = 0; row < 7 ; row++) {
      for (int col = 0; col < 7; col++) {
        font_bw_img[col][row][12] = arr_2d_plus[col][row];
      }
    }

    int arr_2d_minus[7][7] = {{0,0,0,0,0,0,0},
                              {0,0,0,0,0,0,0},
                              {0,0,0,0,0,0,0},
                              {0,1,1,1,1,1,0},
                              {0,0,0,0,0,0,0},
                              {0,0,0,0,0,0,0},
                              {0,0,0,0,0,0,0}};
    for (int row = 0; row < 7 ; row++) {
      for (int col = 0; col < 7; col++) {
        font_bw_img[col][row][13] = arr_2d_minus[col][row];
      }
    }

    int arr_2d_space[7][7] = {{0,0,0,0,0,0,0},
                              {0,0,0,0,0,0,0},
                              {0,0,0,0,0,0,0},
                              {0,0,0,0,0,0,0},
                              {0,0,0,0,0,0,0},
                              {0,0,0,0,0,0,0},
                              {0,0,0,0,0,0,0}};
    for (int row = 0; row < 7 ; row++) {
      for (int col = 0; col < 7; col++) {
        font_bw_img[col][row][14] = arr_2d_space[col][row];
      }
    }

    int arr_2d_null[7][7] =  {{0,0,0,0,0,0,0},
                              {0,0,0,0,0,0,0},
                              {0,0,0,0,0,0,0},
                              {0,0,0,0,0,0,0},
                              {0,0,0,0,0,0,0},
                              {0,0,0,0,0,0,0},
                              {0,0,0,0,0,0,0}};
    for (int row = 0; row < 7 ; row++) {
      for (int col = 0; col < 7; col++) {
        font_bw_img[col][row][15] = arr_2d_null[col][row];
      }
    }

    return err_cod;
  }

  int get_digits( double nm, int (&dgt_num)[12]){

    int err_cod = 0;

    char asc_str[12];

    std::string str;
    std::cout << "nm =" << nm << "\n";
    sprintf (asc_str, "           ");

    sprintf (asc_str, "%g", nm);
    std::cout << "asc_str = <<" << asc_str << ">>\n\n";
    str = asc_str;

    for (int i = 0; i < 12; i++){
      str = asc_str[i];
      dgt_num[i] = asc_str[i] - 48;

      if( dgt_num[i] < 0 or dgt_num[i] > 9 ){
        if( str == "." ){
          dgt_num[i] = 10;
        }else if( str == "e" ){
          dgt_num[i] = 11;
        }else if( str == "+" ){
          dgt_num[i] = 12;
        }else if( str == "-" ){
          dgt_num[i] = 13;
        }else if( str == " " ){
          dgt_num[i] = 14;
        }else if( asc_str[i] == 0 ){
          dgt_num[i] = 15;
        }else{
          std::cout << "\nfound \'" << str << "\' and not converted";
          err_cod = asc_str[i];
          std::cout << "\n char(" << err_cod << ") found\n";
          err_cod = 1;
        }
      }

    //std::cout << " {" << asc_str[i] << "} " << "["<< dgt_num[i] <<"],";
    }

    //std::cout << "\nerr_cod =" << err_cod << "\n";

    return err_cod;
  }


#endif
