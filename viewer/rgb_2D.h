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
#include <string>
#include <cmath>
#include <scitbx/array_family/flex_types.h>

#include <scitbx/array_family/tiny_types.h>

#include <dials/viewer/fonts_2D.h>
#include <dials/viewer/mask_bmp_2D.h>
#include <dials/model/data/mask_code.h>

namespace dials { namespace viewer { namespace boost_python {
  using scitbx::af::flex_double;
  using scitbx::af::flex_int;
  using scitbx::af::flex_grid;

  using scitbx::af::int3;
  using scitbx::af::int6;

  using dials::model::Valid;
  using dials::model::Background;
  using dials::model::BackgroundUsed;
  using dials::model::Foreground;

  /*

  flex_int gen_str_tst(flex_double & data_num) {

    flex_int bmp_dat(flex_grid<>(14, 7, 10),0);
    int npos=data_num.accessor().all()[0];
    int digit_val[15];
    int err_conv = 0;
    double dl_nm;

    for (int pos = 0; pos < npos; pos++) {
      dl_nm = data_num(pos, 0);

      err_conv = get_digits(dl_nm, digit_val);
      if(err_conv == 1){
        std::cout << "\nerror converting\n";
      }

    }

    return bmp_dat;
  }


  flex_int gen_font_img(flex_double & data2d) {

    int ndept=data2d.accessor().all()[0];
    flex_int font_3d_img(flex_grid<>(14, 7, ndept),0);
    int err_conv = 0;
    int font_vol[14][7][16];
    err_conv = get_font_img_array(font_vol);

    //std::cout << "\n ndept =" << ndept << "\n";
    if(err_conv == 0){
      for (int dept = 0; dept < ndept; dept ++){
        for (int col = 0; col < 7; col++){
          for (int row = 0; row < 14; row++){
            //std::cout << "row,col,k =" << row << ", " << col << ", " << dept << "\n";
            font_3d_img(row,col,dept) = font_vol[row][col][dept];
          }
        }
      }
    }

    return font_3d_img;
  }

  */


  class rgb_img
  {

    private:
      int red_byte[255 * 3 + 1];
      int green_byte[255 * 3 + 1];
      int blue_byte[255 * 3 + 1];

      int err_conv;

      int font_vol[14][7][16];
      int mask_vol[85][85][4];

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

        blue_byte[765] = 255;
        green_byte[765] = 255;
        red_byte[765] = 255;

        err_conv = get_font_img_array(font_vol);
        if(err_conv != 0){
          std::cout << "\n ERROR Building fonts internally \n error code ="
                    << err_conv << "\n";
        }
        err_conv = get_mask_img_array(mask_vol);
        if(err_conv != 0){
          std::cout << "\n ERROR Building mask internally \n error code ="
                    << err_conv << "\n";
        }

      }

      flex_int gen_bmp(flex_double & data2d, flex_double & mask2d) {

        int nrow=data2d.accessor().all()[0];
        int ncol=data2d.accessor().all()[1];

        double max = 1, min = -1, loc_cel, dif = 0;

        flex_double scaled_array(flex_grid<>(nrow, ncol),0);

         for (int col = 0; col < ncol; col++) {
          for (int row = 0; row < nrow ; row++) {


            loc_cel = data2d(row, col);
            if(row == 0 and col == 0){
              max = loc_cel;
              min = loc_cel;
            }else{
              if(loc_cel > max){
                max = loc_cel;
              }
              if(loc_cel < min){
                min = loc_cel;
              }
            }
          }
        }

        //std::cout << "\n max = "<< max << "\n";
        //std::cout << "\n min = "<< min << "\n \n";

        dif = max - min;

        for (int col = 0; col < ncol; col++) {
          for (int row = 0; row < nrow; row++) {
            scaled_array(row, col) = 255.0 * 3 * (( data2d(row, col) - min)
                                                          / dif );
          }
        }
        int px_scale = 0;

        if(ncol < 200 and nrow < 200){

          float diagn;
            if(ncol < 20 and nrow < 20){
              px_scale = 85;
              //std::cout << "less than (20 * 20) pixels \n";
            }else{
              diagn = sqrt(ncol * ncol + nrow * nrow);
              px_scale = (1000.0 / diagn);
              //std::cout << "scale = " << px_scale << "\n";
            }
        }else{
            px_scale = 1;
        }

        flex_int bmp_dat(flex_grid<>(nrow * px_scale, ncol * px_scale, 3),0);

        int digit_val[15];
        int pix_row;
        int pix_col;

        for (int col = 0; col < ncol; col++) {
          for (int row = 0; row < nrow; row++) {

            if(px_scale > 1){

              //painting the scaled pixel with the *hot* rowour convention

              for(pix_col = col * px_scale;
                  pix_col < col * px_scale + px_scale;
                  pix_col++){

                for(pix_row = row * px_scale;
                    pix_row < row * px_scale + px_scale;
                    pix_row++){

                  bmp_dat(pix_row, pix_col, 0) = red_byte[int(
                                                 scaled_array(row, col))];
                  bmp_dat(pix_row, pix_col, 1) = green_byte[int(
                                                 scaled_array(row, col))];
                  bmp_dat(pix_row, pix_col, 2) = blue_byte[int(
                                                 scaled_array(row, col))];
                }
              }



              // Painting mask into the scaled pixel
              for(int mask_pix_col = 0, pix_col = col * px_scale;
                  mask_pix_col < px_scale;
                  pix_col++,
                  mask_pix_col++){

                for(int mask_pix_row = 0, pix_row = row * px_scale;
                    mask_pix_row < px_scale;
                    pix_row++,
                    mask_pix_row++){
                  if( mask_vol[mask_pix_row][mask_pix_col][0] == 1 and



                     mask2d(row, col) > 0
                  ){


                    /*
                     *
                      (mask2d(row, col) & Valid == Valid)

                      *
                    or
                      mask_vol[mask_pix_row][mask_pix_col][1] == 1 or
                      mask_vol[mask_pix_row][mask_pix_col][2] == 1 or
                      mask_vol[mask_pix_row][mask_pix_col][3] == 1 ){



        if ((mask[i] & Foreground) == Foreground) {
          if ((mask[i] & Valid) == Valid) {
            sum_p_ += signal[i];
            sum_b_ += background[i];
            n_signal_++;
          } else {
            success_ = false;
          }
        } else if ((mask[i] & bg_code) == bg_code) {
          n_background_++;
        }
                     */

                    bmp_dat(pix_row, pix_col, 0) = 150;
                    bmp_dat(pix_row, pix_col, 1) = 150;
                    bmp_dat(pix_row, pix_col, 2) = 150;
                  }

                }
              }


              // Painting intensity value into the scaled pixel
              err_conv = get_digits(data2d(row, col), digit_val);
              if(err_conv == 0){
                //std::cout << "data2d(row, col) = " << data2d(row, col) << "\n";
                for(int dg_num = 0;
                    dg_num < 12 and digit_val[dg_num] != 15;
                    dg_num++){

                  for(int font_pix_col = 0, pix_col = col * px_scale + dg_num * 7;
                      font_pix_col < 7;
                      pix_col++,
                      font_pix_col++){

                    for(int font_pix_row = 0, pix_row = row * px_scale + 14;
                        font_pix_row < 14;
                        pix_row++,
                        font_pix_row++){

                      if(font_vol[font_pix_row][font_pix_col][digit_val[dg_num]] == 1){
                        if( scaled_array(row, col) < 255 ){
                          bmp_dat(pix_row, pix_col, 0) = 255;
                          bmp_dat(pix_row, pix_col, 1) = 255;
                          bmp_dat(pix_row, pix_col, 2) = 0;
                        }else if( scaled_array(row, col) > 255 * 2 ){
                          bmp_dat(pix_row, pix_col, 0) = 0;
                          bmp_dat(pix_row, pix_col, 1) = 0;
                          bmp_dat(pix_row, pix_col, 2) = 0;
                        }else{
                          bmp_dat(pix_row, pix_col, 0) = 00;
                          bmp_dat(pix_row, pix_col, 1) = 00;
                          bmp_dat(pix_row, pix_col, 2) = 255;
                        }
                      }
                    }
                  }
                }
              }

            }else{
              bmp_dat(row, col, 0) = red_byte[int(scaled_array(row, col))];
              bmp_dat(row, col, 1) = green_byte[int(scaled_array(row, col))];
              bmp_dat(row, col, 2) = blue_byte[int(scaled_array(row, col))];
            }



          }
        }
        //std::cout << "\n here 04 \n";

      return bmp_dat;
      }

  };


}}}

#endif
