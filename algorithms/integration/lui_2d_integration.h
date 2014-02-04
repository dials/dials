/*
 * lui_2d_integration.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: Luis Fuentes-Montero
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_LUI_INTEGRATION_2D_H
#define DIALS_ALGORITHMS_LUI_INTEGRATION_2D_H
#include <stdio.h>
#include <iostream>
#include <scitbx/vec2.h>
#include <cmath>
#include <scitbx/array_family/flex_types.h>

#include <cstdlib>
#include <scitbx/array_family/versa_matrix.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/model/data/shoebox.h>
#include <scitbx/array_family/accessors/mat_grid.h>

namespace dials { namespace algorithms {

  using scitbx::vec2;
  using scitbx::af::flex_int;
  using scitbx::af::flex_double;
  using scitbx::af::flex_grid;
  using dials::model::Foreground;
  using dials::model::Background;
  // using dials::model::Valid;

  vec2<double> raw_2d_cut(
    const af::const_ref< double, af::c_grid<2> > &data2d,
    const af::const_ref< int, af::c_grid<2> > &mask2d,
    const af::const_ref< double, af::c_grid<2> > &background2d) {

    // classic and simple integration summation

    double i_s = 0, i_bg = 0, rho_j = 0;
    double n = 0, m = 0;
    int cont = 0;
    double var_i;
    std::size_t ncol=data2d.accessor()[1];
    std::size_t nrow=data2d.accessor()[0];
    vec2<double> integr_data(0,1);

    // looping thru each pixel
    for (int row = 0; row<nrow;row++) {
      for (int col = 0; col<ncol;col++) {
        if ( mask2d(row,col) & Foreground ) {
          i_s  += data2d(row,col) - background2d(row,col);
          i_bg += background2d(row,col);
          m++;
        } else if(mask2d(row,col) & Background) {
          rho_j += background2d(row,col);
          n++;
        }
        cont++;
      }
    }
    if( i_bg>0 && cont>0 && m>0 && n>0 ){
      // Calculation of variance, by following:
      // A. G. W. Leslie
      // Acta Cryst. (1999). D55, 1696-1702
      // eq # 9
      var_i = i_s + i_bg + (m / n) * ( m / n) * rho_j;
    } else {
      var_i = i_s;
    }

    integr_data[0]=i_s;            // intensity summation
    integr_data[1]=var_i;          // intensity variance
    return integr_data;

  }

  af::versa< double, af::c_grid<2> > add_2d(
    const af::const_ref< double, af::c_grid<2> > &descriptor,
    const af::const_ref< double, af::c_grid<2> > &data2d,
    const af::const_ref< double, af::c_grid<2> > &tmp_total) {
    int ncol_in=data2d.accessor()[1];
    int nrow_in=data2d.accessor()[0];
    int ncol_tot=tmp_total.accessor()[1];
    int nrow_tot=tmp_total.accessor()[0];

    if( ncol_tot < ncol_in or nrow_tot < nrow_in ){
      std::cout << "\n WRONG SISE OF ARRAYS \n";
    }

    af::versa< double, af::c_grid<2> > total(af::c_grid<2>(nrow_tot,ncol_tot),0);

    for (int row = 0; row < nrow_tot; row++) {
      for (int col = 0; col < ncol_tot; col++) {
        total(row,col) = tmp_total(row,col);
      }
    }

    //double total[nrow_tot, ncol_tot];
    //flex_double total(tmp_total);
    double centr_col = descriptor(0,0);
    double centr_row = descriptor(0,1);
    double scale = descriptor(0,2);
    double tot_row, tot_col;
    double x_contrib, y_contrib, x_pix_pos, y_pix_pos;
    int xpos_ex, ypos_ex;
    int tot_row_centr=int(nrow_tot / 2);
    int tot_col_centr=int(ncol_tot / 2);
    // looping thru each pixel
    for (int row = 0; row < nrow_in; row++) {
      for (int col = 0; col < ncol_in; col++) {
        tot_row = row + tot_row_centr - centr_row + 1;
        tot_col = col + tot_col_centr - centr_col + 1;
        if (tot_row > 0 and tot_col > 0 and
          tot_row < nrow_tot and tot_col < ncol_tot) {

          // interpolating and finding contributions of each pixel
          // to the four surrounding  pixels in the other shoe box
          x_pix_pos=tot_col - int(tot_col);
          if (x_pix_pos < 0.5) {
            x_contrib = 0.5 + x_pix_pos;
            xpos_ex = tot_col - 1;
          } else if ( x_pix_pos > 0.5 ) {
            x_contrib = 1.5 - x_pix_pos;
            xpos_ex = tot_col + 1;
          } else {
            x_contrib = 1;
            xpos_ex = tot_col;
          }

          y_pix_pos=tot_row - int(tot_row);
          if (y_pix_pos < 0.5) {
            y_contrib = 0.5 + y_pix_pos;
            ypos_ex = tot_row - 1;
          } else if ( y_pix_pos > 0.5 ) {
            y_contrib = 1.5 - y_pix_pos;
            ypos_ex = tot_row + 1;
          } else {
            y_contrib = 1;
            ypos_ex = tot_row;
          }
          int pos_tot_row = int(tot_row);
          int pos_tot_col = int(tot_col);
          // Adding corresponding contributions to each pixel
          total(pos_tot_row, pos_tot_col) += data2d(row, col) * scale
            * x_contrib * y_contrib;
          if( xpos_ex != tot_col or ypos_ex != tot_row ){
            if( xpos_ex != tot_col ){
              total(pos_tot_row, xpos_ex) += data2d(row,col) * scale
                * (1 - x_contrib) * y_contrib;
            }
            if( ypos_ex != tot_row ){
              total(ypos_ex, pos_tot_col) += data2d(row, col) * scale
                * x_contrib * (1 - y_contrib);
            }
            if( xpos_ex != tot_col and ypos_ex != tot_row ){
              total(ypos_ex, xpos_ex) += data2d(row, col)* scale
                * (1 - x_contrib) * (1 - y_contrib);
            }
          }

        } else {
          std::cout << "\n ERROR Not fitting in the area to be added \n";
          std::cout <<"  tot_row =" <<tot_row << "  tot_col =" << tot_col <<
              "  nrow_tot =" << nrow_tot << "  ncol_tot =" << ncol_tot << "\n";
          std::cout <<" centr_col =" <<  centr_col <<
            " centr_row =" << centr_row << "\n";
        }
      }
    }


    return total;
  }

  double m_linear_scale(int & cnt, double i_mod[], double i_exp[]){
    double m = 0;
    if (cnt > 0){
      double i_wgt[cnt];
      double avg = 0;
      for (int i = 0; i < cnt; i++){
        avg += i_mod[i];
      }
      avg = avg/double(cnt);
      for (int i = 0; i < cnt; i++){
        i_wgt[i] = i_mod[i] / avg;
      }
      avg = 0;
      for (int i = 0; i < cnt; i++){
        avg += i_wgt[i];
      }
      avg = avg/double(cnt);
      double scale = 0, px_scl;
      for (int i = 0; i < cnt; i++){
        px_scl = i_exp[i] / i_mod[i];
        scale += px_scl * i_wgt[i];
      }
      m = scale / double(cnt);
    } else {
      std::cout << "\n missing useful pixels for profile fitting\n";
    }

    return m;
  }

  double w_m_least_squres_1d(int & cnt, double m_in, double i_mod[], double i_exp[]){

    double sum_xy = 0, sum_x_sq = 0, m;
    for (int i = 0; i < cnt; i++){

      sum_xy += i_exp[i];// / m_in;
      sum_x_sq += i_mod[i];// / m_in;
    }
    m = sum_xy / sum_x_sq;
    return m;
  }

  /*
  double m_least_squres_1d(int & cnt, double i_mod[], double i_exp[]){
    // least-squares scaling following the formula:
    // m = ( sum(X(i) * Y(i) ) / sum( X(i)**2) )

    double sum_xy = 0, sum_x_sq = 0, m;
    for (int i = 0; i < cnt; i++){
      sum_xy += i_exp[i] * i_mod[i];
      sum_x_sq += i_mod[i] * i_mod[i];
    }
    m = sum_xy / sum_x_sq;
    return m;
  }
  */


  // Given a 2D shoebox and a 2D profile, fits the profile to find the scale
  vec2<double> fitting_2d(
    const af::const_ref< double, af::c_grid<2> > &descriptor,
    const af::const_ref< double, af::c_grid<2> > &data2d,
    const af::const_ref< double, af::c_grid<2> > &backg2d,
    const af::const_ref< double, af::c_grid<2> > &profile2d) {

    int ncol = profile2d.accessor()[1];
    int nrow = profile2d.accessor()[0];
    int counter = 0;

    af::versa< double, af::c_grid<2> > data2dmov_01(af::c_grid<2>(nrow, ncol),0);
    af::versa< double, af::c_grid<2> > backg2dmov_01(af::c_grid<2>(nrow, ncol),0);

    af::versa< double, af::c_grid<2> > data2dmov(af::c_grid<2>(nrow, ncol),0);
    af::versa< double, af::c_grid<2> > backg2dmov(af::c_grid<2>(nrow, ncol),0);

    //descriptor(0,2) = 1;
    data2dmov = add_2d(descriptor, data2d, data2dmov_01.const_ref());
    backg2dmov = add_2d(descriptor, backg2d, backg2dmov_01.const_ref());

    // Counting how many pixels are useful so far
    for (int row = 0; row < nrow; row++) {
      for (int col = 0; col < ncol; col++) {
        if (data2dmov(row,col) != backg2dmov(row,col) ) {
          counter++ ;
        }
      }
    }

    // Building a set 1D lists with the useful pixels

    double iexpr_lst[counter];
    double imodl_lst[counter];
    double sum = 0, i_var;

    counter = 0;
    for (int row = 0; row < nrow; row++) {
      for (int col = 0; col < ncol; col++) {
        if (data2dmov(row,col) != backg2dmov(row,col)
          and profile2d(row,col) > 0 ) {

          iexpr_lst[counter] = data2dmov(row,col) - backg2dmov(row,col);
          imodl_lst[counter] = profile2d(row,col);// * conv_scale;
          counter++ ;
        }
      }
    }

    // finding the scale needed to fit profile list to experiment list
    double m, diff, df_sqr;


    // this is how the algorithm for fitting must be chosen
    // between linear or least squares

    m = m_linear_scale(counter, imodl_lst, iexpr_lst);
    //m = m_least_squres_1d(counter, imodl_lst, iexpr_lst);
    double tmp_m = m;
    m = w_m_least_squres_1d(counter, tmp_m, imodl_lst, iexpr_lst);


    //measuring R
    sum = 0;
    for (int row = 0; row < nrow; row++) {
      for (int col = 0; col < ncol; col++) {

        diff = profile2d(row, col) * m - data2dmov(row, col);
        df_sqr = diff * diff;
        sum += df_sqr;

      }
    }
    i_var = sum;

    //measuring the volume of the scaled profile
    sum = 0;
    for (int row = 0; row < nrow; row++) {
      for (int col = 0; col < ncol; col++) {
        sum += profile2d(row, col) * m;
      }
    }
    vec2<double> integr_data(0,1);
    integr_data[0] = sum;                   // intensity
    integr_data[1] = i_var;                 // intensity variance

    return integr_data;
  }

  flex_double subtrac_bkg_2d(flex_double data2d, flex_double backg2d) {
    // given a 2D shoebox and a 2D background,
    // it subtract the background from the shoebox
    int ncol = data2d.accessor().all()[1];
    int nrow = data2d.accessor().all()[0];
    flex_double data2dreturn(data2d.accessor(), 0);
    for (int row = 0; row < nrow; row++) {
      for (int col = 0; col < ncol; col++) {
        if (data2d(row,col) > backg2d(row,col) ) {
          data2dreturn(row,col) = data2d(row,col) - backg2d(row,col);
        } else {
          data2dreturn(row,col) = 0.0;
        }
      }
    }
    return data2dreturn;
  }


} }

#endif
