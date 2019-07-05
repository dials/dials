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

#include <dials/viewer/fonts_2D.h>
#include <dials/viewer/mask_bmp_2D.h>
#include <dials/model/data/mask_code.h>

namespace dials { namespace viewer { namespace boost_python {
  using scitbx::af::flex_double;
  using scitbx::af::flex_grid;
  using scitbx::af::flex_int;

  using dials::model::Background;
  using dials::model::BackgroundUsed;
  using dials::model::Foreground;
  using dials::model::Valid;

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

  class rgb_img {
  private:
    int hot_pal_red_byte[255 * 3 + 1];
    int hot_pal_green_byte[255 * 3 + 1];
    int hot_pal_blue_byte[255 * 3 + 1];

    int gray_all_rgb_byte[255 * 3 + 1];

    double max, min;
    int err_conv;

    int font_vol[14][7][16];
    int mask_vol[85][85][4];

  public:
    rgb_img() {
      // Generating grayscale palette
      for (int i = 0; i < 255 * 3; i++) {
        gray_all_rgb_byte[i] = int(i / 3);
      }
      gray_all_rgb_byte[765] = 255;

      // Generating hot colour palette
      for (int i = 0; i < 255; i++) {
        hot_pal_red_byte[i] = i;
        hot_pal_green_byte[i + 255] = i;
        hot_pal_blue_byte[i + 255 * 2] = i;
      }

      for (int i = 255; i < 255 * 3; i++) {
        hot_pal_red_byte[i] = 255;
      }

      for (int i = 0; i < 255; i++) {
        hot_pal_green_byte[i] = 0;
      }

      for (int i = 255 * 2; i < 255 * 3; i++) {
        hot_pal_green_byte[i] = 255;
      }

      for (int i = 0; i < 255 * 2; i++) {
        hot_pal_blue_byte[i] = 0;
      }

      hot_pal_blue_byte[765] = 255;
      hot_pal_green_byte[765] = 255;
      hot_pal_red_byte[765] = 255;

      max = -1;
      min = -1;

      err_conv = get_font_img_array(font_vol);
      if (err_conv != 0) {
        std::cout << "\n ERROR Building fonts internally \n error code =" << err_conv
                  << "\n";
      }
      err_conv = get_mask_img_array(mask_vol);
      if (err_conv != 0) {
        std::cout << "\n ERROR Building mask internally \n error code =" << err_conv
                  << "\n";
      }
    }

    int set_min_max(double new_min, double new_max) {
      if (new_min < new_max) {
        min = new_min;
        max = new_max;
      } else {
        min = new_min;
        max = min + 1;
      }

      // std::cout << "\n min(new), max(new) =" << min << ", " << max << "\n";
      return 0;
    }

    flex_int gen_bmp(flex_double& data2d,
                     flex_double& mask2d,
                     bool show_nums,
                     int palette_num) {
      // debugging palette number passed from Python
      /*
      std::cout << "\n show_nums =" << show_nums << "\n";
      std::cout << "\n palette_num =" << palette_num << "\n";
      */

      // Python code that enumerates palette
      /*
          print "palette =", palette
          if( palette == "black2white" ):
            palette_num = 1
          elif( palette == "white2black" ):
            palette_num = 2
          elif( palette == "hot ascend" ):
            palette_num = 3
          else: # assuming "hot descend"
            palette_num = 4
      */

      int nrow = data2d.accessor().all()[0];
      int ncol = data2d.accessor().all()[1];
      int col, row;
      double loc_cel, dif = 0;
      int loc_cel_int;

      // flex_double scaled_array(flex_grid<>(nrow, ncol),0);

      if (max == -1 && min == -1) {
        for (row = 0; row < nrow; row++) {
          for (col = 0; col < ncol; col++) {
            loc_cel = data2d(row, col);
            if (row == 0 && col == 0) {
              max = loc_cel;
              min = loc_cel;
            } else {
              if (loc_cel > max) {
                max = loc_cel;
              }
              if (loc_cel < min) {
                min = loc_cel;
              }
            }
          }
        }
      }

      dif = max - min;

      int px_scale = 0;

      if (ncol < 200 && nrow < 200) {
        px_scale = 85;
        // std::cout << "less than (200 * 200) pixels \n";
      } else {
        px_scale = 1;
      }
      // std::cout << "\n px_scale = " << px_scale << "\n";

      flex_int bmp_dat(flex_grid<>(nrow * px_scale, ncol * px_scale, 3), 0);

      int digit_val[15];
      int pix_row, pix_col;

      int mask_pix_col, mask_pix_row;
      double scaled_pixel;

      // std::cout << "\n ncol =" << ncol << " \n";
      // std::cout << "\n nrow =" << nrow << " \n";
      // std::cout << "\n palette_num =" << palette_num << "\n";

      for (row = 0; row < nrow; row++) {
        for (col = 0; col < ncol; col++) {
          loc_cel = data2d(row, col);
          if (loc_cel > max) {
            loc_cel = max;
          }
          if (loc_cel < min) {
            loc_cel = min;
          }
          scaled_pixel = 255.0 * 3 * ((loc_cel - min) / dif);

          loc_cel_int = int(mask2d(row, col));

          if (px_scale > 1) {
            // painting the scaled pixel with the *hot* color convention
            for (pix_col = col * px_scale; pix_col < col * px_scale + px_scale;
                 pix_col++) {
              for (pix_row = row * px_scale; pix_row < row * px_scale + px_scale;
                   pix_row++) {
                if (palette_num == 1) {
                  bmp_dat(pix_row, pix_col, 0) = gray_all_rgb_byte[int(scaled_pixel)];
                  bmp_dat(pix_row, pix_col, 1) = gray_all_rgb_byte[int(scaled_pixel)];
                  bmp_dat(pix_row, pix_col, 2) = gray_all_rgb_byte[int(scaled_pixel)];
                } else if (palette_num == 2) {
                  bmp_dat(pix_row, pix_col, 0) =
                    gray_all_rgb_byte[int(765 - scaled_pixel)];
                  bmp_dat(pix_row, pix_col, 1) =
                    gray_all_rgb_byte[int(765 - scaled_pixel)];
                  bmp_dat(pix_row, pix_col, 2) =
                    gray_all_rgb_byte[int(765 - scaled_pixel)];
                } else if (palette_num == 3) {
                  bmp_dat(pix_row, pix_col, 0) = hot_pal_red_byte[int(scaled_pixel)];
                  bmp_dat(pix_row, pix_col, 1) = hot_pal_green_byte[int(scaled_pixel)];
                  bmp_dat(pix_row, pix_col, 2) = hot_pal_blue_byte[int(scaled_pixel)];
                } else {
                  bmp_dat(pix_row, pix_col, 0) =
                    hot_pal_red_byte[int(765 - scaled_pixel)];
                  bmp_dat(pix_row, pix_col, 1) =
                    hot_pal_green_byte[int(765 - scaled_pixel)];
                  bmp_dat(pix_row, pix_col, 2) =
                    hot_pal_blue_byte[int(765 - scaled_pixel)];
                }
              }
            }

            // Painting mask into the scaled pixel
            for (mask_pix_col = 0, pix_col = col * px_scale; mask_pix_col < px_scale;
                 pix_col++, mask_pix_col++) {
              for (mask_pix_row = 0, pix_row = row * px_scale; mask_pix_row < px_scale;
                   pix_row++, mask_pix_row++) {
                if ((mask_vol[mask_pix_row][mask_pix_col][0] == 1
                     && ((loc_cel_int & Valid) == Valid))
                    || (mask_vol[mask_pix_row][mask_pix_col][1] == 1
                        && ((loc_cel_int & Foreground) == Foreground))

                    || (mask_vol[mask_pix_row][mask_pix_col][2] == 1
                        && ((loc_cel_int & BackgroundUsed) == BackgroundUsed))
                    || (mask_vol[mask_pix_row][mask_pix_col][3] == 1
                        && ((loc_cel_int & Background) == Background))) {
                  if (palette_num == 1 || palette_num == 2) {
                    bmp_dat(pix_row, pix_col, 0) = 250;
                    bmp_dat(pix_row, pix_col, 1) = 50;
                    bmp_dat(pix_row, pix_col, 2) = 50;
                  } else {
                    bmp_dat(pix_row, pix_col, 0) = 150;
                    bmp_dat(pix_row, pix_col, 1) = 150;
                    bmp_dat(pix_row, pix_col, 2) = 150;
                  }
                }
              }
            }

            // Painting intensity value into the scaled pixel

            if (show_nums == true) {
              err_conv = get_digits(data2d(row, col), digit_val);
              if (err_conv == 0) {
                // std::cout << "data2d(row, col) = " << data2d(row, col) << "\n";
                for (int dg_num = 0; dg_num < 12 && digit_val[dg_num] != 15; dg_num++) {
                  for (int font_pix_col = 0, pix_col = col * px_scale + dg_num * 7;
                       font_pix_col < 7;
                       pix_col++, font_pix_col++) {
                    for (int font_pix_row = 0, pix_row = row * px_scale + 14;
                         font_pix_row < 14;
                         pix_row++, font_pix_row++) {
                      if (font_vol[font_pix_row][font_pix_col][digit_val[dg_num]]
                          == 1) {
                        if (palette_num == 1 || palette_num == 3) {
                          if (scaled_pixel < 255) {
                            bmp_dat(pix_row, pix_col, 0) = 255;
                            bmp_dat(pix_row, pix_col, 1) = 255;
                            bmp_dat(pix_row, pix_col, 2) = 0;
                          } else if (scaled_pixel > 255 * 2) {
                            bmp_dat(pix_row, pix_col, 0) = 0;
                            bmp_dat(pix_row, pix_col, 1) = 0;
                            bmp_dat(pix_row, pix_col, 2) = 0;
                          } else {
                            bmp_dat(pix_row, pix_col, 0) = 00;
                            bmp_dat(pix_row, pix_col, 1) = 00;
                            bmp_dat(pix_row, pix_col, 2) = 255;
                          }
                        } else {
                          if (scaled_pixel < 255) {
                            bmp_dat(pix_row, pix_col, 0) = 0;
                            bmp_dat(pix_row, pix_col, 1) = 0;
                            bmp_dat(pix_row, pix_col, 2) = 0;
                          } else if (scaled_pixel > 255 * 2) {
                            bmp_dat(pix_row, pix_col, 0) = 255;
                            bmp_dat(pix_row, pix_col, 1) = 255;
                            bmp_dat(pix_row, pix_col, 2) = 0;
                          } else {
                            bmp_dat(pix_row, pix_col, 0) = 00;
                            bmp_dat(pix_row, pix_col, 1) = 00;
                            bmp_dat(pix_row, pix_col, 2) = 255;
                          }
                        }
                      }
                    }
                  }
                }
              }
            }

          } else {
            if (palette_num == 1) {
              bmp_dat(row, col, 0) = gray_all_rgb_byte[int(scaled_pixel)];
              bmp_dat(row, col, 1) = gray_all_rgb_byte[int(scaled_pixel)];
              bmp_dat(row, col, 2) = gray_all_rgb_byte[int(scaled_pixel)];
            } else if (palette_num == 2) {
              bmp_dat(row, col, 0) = gray_all_rgb_byte[int(765 - scaled_pixel)];
              bmp_dat(row, col, 1) = gray_all_rgb_byte[int(765 - scaled_pixel)];
              bmp_dat(row, col, 2) = gray_all_rgb_byte[int(765 - scaled_pixel)];
            } else if (palette_num == 3) {
              bmp_dat(row, col, 0) = hot_pal_red_byte[int(scaled_pixel)];
              bmp_dat(row, col, 1) = hot_pal_green_byte[int(scaled_pixel)];
              bmp_dat(row, col, 2) = hot_pal_blue_byte[int(scaled_pixel)];
            } else {
              bmp_dat(row, col, 0) = hot_pal_red_byte[int(765 - scaled_pixel)];
              bmp_dat(row, col, 1) = hot_pal_green_byte[int(765 - scaled_pixel)];
              bmp_dat(row, col, 2) = hot_pal_blue_byte[int(765 - scaled_pixel)];
            }
          }
        }
      }
      // std::cout << "\n building BMP \n";

      return bmp_dat;
    }
  };

}}}  // namespace dials::viewer::boost_python

#endif
