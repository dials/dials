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
  using dials::model::Valid;
  using dials::model::Background;
  using std::sqrt;
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

    integr_data[0] = i_s;            // intensity summation
    integr_data[1] = var_i;          // intensity variance
    return integr_data;

  }

  double sigma_2d(
    const float intensity,
    const af::const_ref< int, af::c_grid<2> > &mask2d,
    const af::const_ref< double, af::c_grid<2> > &background2d) {

    // classic and simple integration summation

    double itst_bkg = 0, rho_j = 0;
    double n = 0, m = 0;
    int cont = 0;
    double var_i;
    std::size_t ncol=background2d.accessor()[1];
    std::size_t nrow=background2d.accessor()[0];
    //vec2<double> integr_data(0,1);
    double integr_data;

    // looping thru each pixel
    for (int row = 0; row<nrow;row++) {
      for (int col = 0; col<ncol;col++) {
        if ( mask2d(row,col) & Foreground ) {

          itst_bkg += background2d(row,col);
          m++;
        } else if(mask2d(row,col) & Background) {
          rho_j += background2d(row,col);
          n++;
        }
        cont++;
      }
    }
    if( itst_bkg>0 && cont>0 && m>0 && n>0 ){
      // Calculation of variance, by following:
      // A. G. W. Leslie
      // Acta Cryst. (1999). D55, 1696-1702
      // eq # 9
      var_i = intensity + itst_bkg + (m / n) * ( m / n) * rho_j;
    } else {
      var_i = intensity;
    }

    integr_data = var_i;          // intensity variance

    return integr_data;

  }
  af::versa< int, af::c_grid<2> > mask_2d_interpolate(
    const af::const_ref< double, af::c_grid<2> > &descriptor,
    const af::const_ref< int, af::c_grid<2> > &mask2d_in,
    const af::const_ref< int, af::c_grid<2> > &mask2d_tmp_total) {

    int nrow_tot = mask2d_tmp_total.accessor()[0];
    int ncol_tot = mask2d_tmp_total.accessor()[1];

    int nrow_in = mask2d_in.accessor()[0];
    int ncol_in = mask2d_in.accessor()[1];

    double tot_row, tot_col;
    //double x_pix_pos, y_pix_pos;
    //int xpos_ex, ypos_ex;
    int tot_row_centr = int(nrow_tot / 2);
    int tot_col_centr = int(ncol_tot / 2);

    double centr_col = descriptor(0,0);
    double centr_row = descriptor(0,1);
    af::versa< int, af::c_grid<2> > mask2d_out
                   (af::c_grid<2>(nrow_tot, ncol_tot), 0);

    for (int row = 0; row < nrow_in; row++) {
      for (int col = 0; col < ncol_in; col++) {
        tot_row = row + tot_row_centr - centr_row + 1;
        tot_col = col + tot_col_centr - centr_col + 1;
        mask2d_out(tot_row, tot_col) = mask2d_in(row, col);
      }
    }

    for (int row = 0; row < nrow_in; row++) {
      for (int col = 0; col < ncol_in; col++) {

        if (mask2d_in(row, col) <= 0){

          tot_row = row + tot_row_centr - centr_row + 0.5;
          tot_col = col + tot_col_centr - centr_col + 0.5;
          mask2d_out(tot_row, tot_col) = mask2d_in(row, col);

          tot_row = row + tot_row_centr - centr_row + 1.5;
          tot_col = col + tot_col_centr - centr_col + 0.5;
          mask2d_out(tot_row, tot_col) = mask2d_in(row, col);

          tot_row = row + tot_row_centr - centr_row + 0.5;
          tot_col = col + tot_col_centr - centr_col + 1.5;
          mask2d_out(tot_row, tot_col) = mask2d_in(row, col);

          tot_row = row + tot_row_centr - centr_row + 1.5;
          tot_col = col + tot_col_centr - centr_col + 1.5;
          mask2d_out(tot_row, tot_col) = mask2d_in(row, col);

        }
      }
    }

    /*
    std::cout << "\n" << "desc(0) ="  << descriptor(0,0) << "\n"
              << "\n" << "desc(1) ="  << descriptor(0,1) << "\n"
              << "\n" << "desc(2) ="  << descriptor(0,2) << "\n";
    */

    return mask2d_out;
  }


  // function used repeatedly for adding and interpolation
  // of 2D images with a reflection on it
  //
  // is possible to position and scale the reflection to be added
  // By managing the values in the descriptor variable

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

          std::cout << " ncol_in =" << ncol_in << "  nrow_in =" << nrow_in << "\n";


          std::cout <<"  tot_row =" <<tot_row << "  tot_col =" << tot_col <<
              "\n  nrow_tot =" << nrow_tot << "  ncol_tot =" << ncol_tot << "\n";

          std::cout << " centr_row =" << centr_row <<
            " centr_col =" <<  centr_col << "\n";
        }
      }
    }


    return total;
  }





  af::versa< double, af::c_grid<2> > simple_2d_add(
    const af::const_ref< double, af::c_grid<2> > &in_data2d_one,
    const af::const_ref< double, af::c_grid<2> > &in_data2d_two) {
    int ncol_one=in_data2d_one.accessor()[1];
    int nrow_one=in_data2d_one.accessor()[0];
    int ncol_two=in_data2d_two.accessor()[1];
    int nrow_two=in_data2d_two.accessor()[0];

    if( ncol_two != ncol_one or nrow_two != nrow_one ){
      std::cout << "\n WRONG SISE OF ARRAYS \n";
      std::cout << "\n ncol_one=" << ncol_one << "\n";
      std::cout << "\n ncol_two=" << ncol_two << "\n";
      std::cout << "\n nrow_one=" << nrow_one << "\n";
      std::cout << "\n nrow_two=" << nrow_two << "\n";
    }

    af::versa< double, af::c_grid<2> > total(af::c_grid<2>(nrow_two,ncol_two),0);

    for (int row = 0; row < nrow_two; row++) {
      for (int col = 0; col < ncol_two; col++) {
        total(row,col) =in_data2d_one(row, col) + in_data2d_two(row, col);
      }
    }

    return total;
  }

  af::versa< int, af::c_grid<2> > mask_add_2d(
    const af::const_ref< int, af::c_grid<2> > &mask2d_one,
    const af::const_ref< int, af::c_grid<2> > &mask2d_two) {

    int ncol=mask2d_one.accessor()[1];
    int nrow=mask2d_one.accessor()[0];
    int ncol_tst=mask2d_two.accessor()[1];
    int nrow_tst=mask2d_two.accessor()[0];

    if( ncol_tst != ncol or nrow_tst != nrow ){
      std::cout << "\n WRONG SISE OF ARRAYS \n";
      std::cout << "\n ncol=" << ncol << "\n";
      std::cout << "\n ncol_tst=" << ncol_tst << "\n";
      std::cout << "\n nrow=" << nrow << "\n";
      std::cout << "\n nrow_tst=" << nrow_tst << "\n";
    }

    af::versa< int, af::c_grid<2> > final_mask(af::c_grid<2>(nrow, ncol),0);

    for (int row = 0; row < nrow; row++) {
      for (int col = 0; col < ncol; col++) {

        //  [ || ] = .or.  returning a boolean
        //  [ && ] = .and. returning a boolean

        //  [ | ] = .or.   as a mathematical operation  bit per bit
        //  [ & ] = .and.  as a mathematical operation  bit per bit

        // [ ! ] = .not. operator

        if ( (mask2d_one(row,col) & Background)
         ||  (mask2d_two(row,col) & Background) )  {
          //final_mask(row, col) = final_mask(row, col) | Background;
          final_mask(row, col) = Background | Valid;
        }
      }
    }

    for (int row = 0; row < nrow; row++) {
      for (int col = 0; col < ncol; col++) {
        if ( (mask2d_one(row,col) & Foreground)
         ||  (mask2d_two(row,col) & Foreground) )  {
          //final_mask(row, col) = final_mask(row, col) | Foreground;
          final_mask(row, col) = Foreground | Valid;

        if ( ( !(mask2d_one(row,col) & Valid ) )
         ||  ( !(mask2d_two(row,col) & Valid ) ) ) {
          //final_mask(row, col) &= ~Valid;
          final_mask(row, col) = 0;
        }

        }
      }
    }


    return final_mask;
  }



  // 1D weighted least squares for partially recorded reflections
  double w_least_squares_1d(int & cnt, double i_mod[],
                            double i_exp[], double w[]){

    double sum_xy = 0, sum_x_sq = 0, m;
    for (int i = 0; i < cnt; i++){
      sum_xy += i_mod[i] * i_exp[i] * w[i];
      sum_x_sq += i_mod[i] * i_mod[i] * w[i];
    }
    m = sum_xy / sum_x_sq;
    return m;
  }

  // Given a 2D shoebox and a 2D profile, fits the profile to find the scale
  // used in partially recorded reflections where should not
  // be refined the background plane
  vec2<double> fitting_2d_partials(
    const af::const_ref< double, af::c_grid<2> > &descriptor,
    //const af::const_ref< double, af::c_grid<2> > &data2d,
    const af::const_ref< double, af::c_grid<2> > &data2dmov,
    //const af::const_ref< double, af::c_grid<2> > &backg2d,
    const af::const_ref< double, af::c_grid<2> > &backg2dmov,
    const af::const_ref< double, af::c_grid<2> > &profile2d,
    const af::const_ref< int, af::c_grid<2> > &interpolation_mask2d,
    double sum_its) {

    int ncol = profile2d.accessor()[1];
    int nrow = profile2d.accessor()[0];
    int counter = 0;

    //af::versa< double, af::c_grid<2> > data2dmov_01(af::c_grid<2>(nrow, ncol),0);
    //af::versa< double, af::c_grid<2> > backg2dmov_01(af::c_grid<2>(nrow, ncol),0);

    //af::versa< double, af::c_grid<2> > data2dmov(af::c_grid<2>(nrow, ncol),0);
    //af::versa< double, af::c_grid<2> > backg2dmov(af::c_grid<2>(nrow, ncol),0);

    //data2dmov = add_2d(descriptor, data2d, data2dmov_01.const_ref());
    //backg2dmov = add_2d(descriptor, backg2d, backg2dmov_01.const_ref());

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
    double w_lst[counter];
    double sum = 0, i_var;
    double predicted_i;

    counter = 0;
    for (int row = 0; row < nrow; row++) {
      for (int col = 0; col < ncol; col++) {

        /*
        if (data2dmov(row,col) != backg2dmov(row,col)
          and profile2d(row,col) > 0 and data2dmov(row,col) > 0
          and backg2dmov(row,col) > 0) {
        */


        if (data2dmov(row,col) != backg2dmov(row,col)
           //and ( interpolation_mask2d(row,col) & Foreground )
           and profile2d(row,col) > 0 and data2dmov(row,col) > 0
           and backg2dmov(row,col) > 0 ) {

          iexpr_lst[counter] = data2dmov(row,col) - backg2dmov(row,col);
          imodl_lst[counter] = profile2d(row,col);// * conv_scale;
          predicted_i = backg2dmov(row,col) + imodl_lst[counter] * sum_its;
          if (predicted_i <= 0){
            predicted_i = 0.000001;
            //std::cout << "\n gone below zero\n";
          }
          w_lst[counter] = 1.0 /
                        (predicted_i);
          counter++;
        }
      }
    }


    // finding the scale needed to fit profile list to experiment list
    double m, diff, df_sqr;

    m = w_least_squares_1d(counter, imodl_lst, iexpr_lst, w_lst);

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

  // Given a 2D shoebox and a 2D profile, builds the matrix to solve
  vec2<double> fitting_2d_multile_var_build_mat(
    const af::const_ref< double, af::c_grid<2> > &descriptor,
    const af::const_ref< double, af::c_grid<2> > &data2d,
    const af::const_ref< double, af::c_grid<2> > &backg2d,
    const af::const_ref< double, af::c_grid<2> > &profile2d,
    double tmp_scale,
    af::ref< double, af::c_grid<2> > mat_a,
    af::ref< double, af::c_grid<2> > vec_b){

    int ncol = profile2d.accessor()[1];
    int nrow = profile2d.accessor()[0];

    af::versa< double, af::c_grid<2> >data2dmov_01(af::c_grid<2>(nrow, ncol),0);
    af::versa< double, af::c_grid<2> >backg2dmov_01(af::c_grid<2>(nrow, ncol),0);

    af::versa< double, af::c_grid<2> > data2dmov(af::c_grid<2>(nrow, ncol),0);
    af::versa< double, af::c_grid<2> > backg2dmov(af::c_grid<2>(nrow, ncol),0);

    //descriptor(0,2) = 1;
    data2dmov = add_2d(descriptor, data2d, data2dmov_01.const_ref());
    backg2dmov = add_2d(descriptor, backg2d, backg2dmov_01.const_ref());

    //elements of the matrix
    double sum_pr_sqr = 0.0;
    double sum_p_pr   = 0.0;
    double sum_q_pr   = 0.0;
    double sum_pr     = 0.0;
    double sum_p_sqr  = 0.0;
    double sum_p_q    = 0.0;
    double sum_p      = 0.0;
    double sum_q_sqr  = 0.0;
    double sum_q      = 0.0;

    double sum_one    = 0.0;

    double sum_pr_ro  = 0.0;
    double sum_p_ro   = 0.0;
    double sum_q_ro   = 0.0;
    double sum_ro     = 0.0;
    double p, q , pr, ro, pd_it, w;

    //looping trough lists and building the elements of the matrix to solve
    for (int row = 0; row<nrow; row++) {
      for (int col = 0; col<ncol; col++) {
        //                          fix me     (the mask should be considered)
        //if (( mask2d(row, col) & Background) or
        //                          fix me     (the mask should be considered)
        //    ( mask2d(row, col) & Foreground) ) {
          p = col + 0.5;
          q = row + 0.5;
          pr = profile2d(row, col);
          ro = data2dmov(row, col);

          pd_it = backg2dmov(row, col);
          if(pr > 0){
            pd_it += pr * tmp_scale;
          }

          if(pd_it < 1.0){
            pd_it = 1.0;
          }

          w = 1.0 / pd_it;

          sum_pr_sqr += w * pr * pr;
          sum_p_pr   += w * p * pr;
          sum_q_pr   += w * q * pr;
          sum_pr     += w * pr;
          sum_p_sqr  += w * p * p;
          sum_p_q    += w * p * q;
          sum_p      += w * p;
          sum_q_sqr  += w * q * q;
          sum_q      += w * q;
          sum_one    += w;

          sum_pr_ro  += w * pr * ro;
          sum_p_ro   += w * p * ro;
          sum_q_ro   += w * q * ro;
          sum_ro     += w * ro;
        //}  // fix me (the mask should be considered)
      }
    }

    mat_a(0,0) = sum_pr_sqr;
    mat_a(0,1) = sum_p_pr;
    mat_a(0,2) = sum_q_pr;
    mat_a(0,3) = sum_pr;

    mat_a(1,0) = sum_p_pr;
    mat_a(1,1) = sum_p_sqr;
    mat_a(1,2) = sum_p_q;
    mat_a(1,3) = sum_p;

    mat_a(2,0) = sum_q_pr;
    mat_a(2,1) = sum_p_q;
    mat_a(2,2) = sum_q_sqr;
    mat_a(2,3) = sum_q;

    mat_a(3,0) = sum_pr;
    mat_a(3,1) = sum_p;
    mat_a(3,2) = sum_q;
    mat_a(3,3) = sum_one;

    vec_b(0,0) = sum_pr_ro;
    vec_b(0,1) = sum_p_ro;
    vec_b(0,2) = sum_q_ro;
    vec_b(0,3) = sum_ro;

    int ok = 0;
  return ok;
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
