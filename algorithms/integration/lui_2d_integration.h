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
#include <scitbx/array_family/flex_types.h>
#include <cmath>
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

    
  flex_double add_2d(flex_double descriptor, flex_double data2d, 
                     flex_double tmp_total) {
    
    // given two shoe boxes, adds the second to the first one,
    // but each pixel interpolates by following the
    // position given in the descriptor
    //
    // This subroutine helps building a profile

    flex_double total(tmp_total);
    int ncol_in = data2d.accessor().all()[1];
    int nrow_in = data2d.accessor().all()[0];
    int ncol_tot = tmp_total.accessor().all()[1];
    int nrow_tot = tmp_total.accessor().all()[0];
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
        if (tot_row >= 0 and tot_col >= 0 and tot_row < nrow_tot and tot_col < ncol_tot) {
          
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
          
          // Adding corresponding contributions to each pixel
          total(tot_row, tot_col) = tmp_total(tot_row, tot_col) + data2d(row, col) * scale * x_contrib * y_contrib;
          if( xpos_ex != tot_col or ypos_ex != tot_row ){
            if( xpos_ex != tot_col ){
              total(tot_row, xpos_ex)=tmp_total(tot_row, xpos_ex) + data2d(row,col) * scale * (1 - x_contrib) * y_contrib;
            }
            if( ypos_ex != tot_row ){
              total(ypos_ex, tot_col)=tmp_total(ypos_ex, tot_col) + data2d(row, col) * scale * x_contrib * (1 - y_contrib);
            }
            if( xpos_ex != tot_col and ypos_ex != tot_row ){
              total(ypos_ex, xpos_ex)=tmp_total(ypos_ex, xpos_ex) + data2d(row, col)* scale * (1 - x_contrib) * (1 - y_contrib);
            }
          }
    
        } else {
          std::cout << "\n ERROR Not fitting in the area to be added \n";
          std::cout <<"  tot_row =" <<tot_row << "  tot_col =" << tot_col <<
              "  nrow_tot =" << nrow_tot << "  ncol_tot =" << ncol_tot << "\n";
          std::cout <<" centr_col =" <<  centr_col << " centr_row =" << centr_row << "\n";
        }
      }
    }
    
    return total;
  }

  

  vec2<double> fitting_2d(flex_double descriptor, flex_double data2d,
                          flex_double backg2d, flex_double profile2d) {
    
    // Given a 2D shoebox and a 2D profile, fits the profile to find the scale

    int ncol = profile2d.accessor().all()[1];
    int nrow = profile2d.accessor().all()[0];
    int counter = 0;
    vec2<double> integr_data(0,1);
    flex_double data2dmov(profile2d.accessor(), 0);
    flex_double backg2dmov(profile2d.accessor(), 0);
    descriptor(0,2) = 1;
    data2dmov = add_2d(descriptor, data2d, data2dmov);
    backg2dmov = add_2d(descriptor, backg2d, backg2dmov);
    
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
    double iback_lst[counter];
    double modl_scal_lst[counter];
    double scale = 0, sum = 0, i_var;
    double avg_i_scale, diff, df_sqr;
    counter = 0;
    for (int row = 0; row < nrow; row++) {
      for (int col = 0; col < ncol; col++) {
        if (data2dmov(row,col) != backg2dmov(row,col) and profile2d(row,col) > 0 ) {
          iexpr_lst[counter] = data2dmov(row,col) - backg2dmov(row,col);
          imodl_lst[counter] = profile2d(row,col);// * conv_scale;
          iback_lst[counter] = backg2dmov(row,col);
          counter++ ;
        }
      }
    }
    
    // finding the scale needed to fit profile list to experiment list
    for (int i = 0; i < counter; i++){
      scale = iexpr_lst[i] / imodl_lst[i];
      sum += scale * imodl_lst[i];
    }
    
    avg_i_scale = sum;
    for (int i = 0; i < counter; i++){
      modl_scal_lst[i] =imodl_lst[i] * avg_i_scale;
    }
    
    sum = 0;
    for (int i = 0; i < counter; i++){
      diff = (modl_scal_lst[i] - iexpr_lst[i]);
      df_sqr = diff * diff;
      sum += df_sqr;
    }
    
    i_var = sum;
    
    sum = 0;
      for (int i = 0; i < counter; i++){
        df_sqr = iback_lst[i] * iback_lst[i];
        sum += df_sqr;
      }
    
    integr_data[0] = avg_i_scale;                   // intensity
    integr_data[1] = i_var;                         // intensity variance
    
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
   
   /*  
    * 
   // this piece of code chould be analized with richard help from Richard
   af::versa< double, af::c_grid<2> > add_2d(
     const af::const_ref< double, af::c_grid<2> > &descriptor,
     const af::const_ref< double, af::c_grid<2> > &data2d,
     const af::const_ref< double, af::c_grid<2> > &tmp_total) {
     af::versa< double, af::c_grid<2> > total(tmp_total.accessor(),0);
     int ncol_in=data2d.accessor()[1];
     int nrow_in=data2d.accessor()[0];
     int ncol_tot=tmp_total.accessor()[1];
     int nrow_tot=tmp_total.accessor()[0];
     //std::cout << "\n ncol_tot =" << ncol_tot << "\n nrow_tot =" << nrow_tot << 
     //"\n total.accessor()[0] =" << total.accessor()[0] << 
     //"\n total.accessor()[1] =" << total.accessor()[1];
     //std::copy(tmp_total.begin(), tmp_total.end(), total.begin());
     for (int row = 0; row < nrow_tot; row++) {
       for (int col = 0; col < ncol_tot; col++) {
         //std::cout << "\n row =" << row << "\n col =" << col << "\n";
         total(row, col) = tmp_total(row, col);
       }
     }    
   
     //printing & testing
     std::cout << "\n <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< printing & testing ";
     for (int row = 0; row < nrow_tot; row++) {
       if (row==0){
         std::cout << "\n  [ [ ";
       } else {
         std::cout << "\n    [ ";
       }
       for (int col = 0; col < ncol_tot;col++) {
         //[ 6 ] = minimum width (no maximum given)
         //[ 2 ] = precision after the period
         printf("%6.3f ", total(row,col));
         //std::cout << int(matx2d[row][col]) << " ,   ";
       }
       //fflush(stdout);
       std::cout << "  ]";
     }
     std::cout << " ] \n";
     
     
     for (int row = 0; row < nrow_tot; row++) {
       if (row==0){
         std::cout << "\n  [ [ ";
       } else {
         std::cout << "\n    [ ";
       }
       for (int col = 0; col < ncol_tot;col++) {
         //[ 6 ] = minimum width (no maximum given)
         //[ 2 ] = precision after the period
         printf("%6.3f ", tmp_total(row,col));
         //std::cout << int(matx2d[row][col]) << " ,   ";
       }
       //fflush(stdout);
       std::cout << "  ]";
     }
     std::cout << " ] \n";
     
     
     for (int row = 0; row < nrow_in; row++) {
       if (row==0){
         std::cout << "\n  [ [ ";
       } else {
         std::cout << "\n    [ ";
       }
       for (int col = 0; col < ncol_in;col++) {
         //[ 6 ] = minimum width (no maximum given)
         //[ 2 ] = precision after the period
         printf("%6.3f ", data2d(row,col));
         //std::cout << int(matx2d[row][col]) << " ,   ";
       }
       //fflush(stdout);
       std::cout << "  ]";
     }
     std::cout << " ] \n";
     
     std::cout << "\n printing & testing >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ";
     
     */  
     
 
} }

#endif
