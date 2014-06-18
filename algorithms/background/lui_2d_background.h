#ifndef DIALS_ALGORITHMS_LUI_2D_BACKGROUND_H
#define DIALS_ALGORITHMS_LUI_2D_BACKGROUND_H
#include <iostream>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/model/data/shoebox.h>

namespace dials { namespace algorithms {
  using dials::model::Foreground;
  using dials::model::Background;
  af::versa< double, af::c_grid<2> > flat_background_flex_2d(
    const af::const_ref< double, af::c_grid<2> > &data2d,
    const af::const_ref< int, af::c_grid<2> > &mask2d) {
        std::size_t ncol=data2d.accessor()[1];
        std::size_t nrow=data2d.accessor()[0];
        af::versa< double, af::c_grid<2> > background2d(data2d.accessor(),0);
        double cont = 0.0, tot_bkgr = 0.0;
        double avg_bkgr = 0;
        for (int row = 0; row < nrow; row++) {
          for (int col = 0; col < ncol; col++) {
            // remember to make sure that a pixel masked as Background
            // cannot be a Foreground pixle (ask James)
            if ( mask2d(row, col) & Background ){
              tot_bkgr += data2d(row, col);
              cont++;
            }
            if ( tot_bkgr > 0 and cont > 0) {
              avg_bkgr = tot_bkgr / cont;
            }
          }
        }
        for (int row = 0; row < nrow; row++) {
          for (int col = 0; col < ncol; col++) {
            if ( mask2d(row, col) & Foreground ){
              background2d(row, col) = avg_bkgr;
            } else {
              background2d(row, col) = data2d(row, col);
            }
          }
        }
    return background2d;
  }


  class local_avg{

    public:
      int max_row, max_col;
      double count, total;

      local_avg(int nm_row, int nm_col){
        max_row = nm_row;
        max_col = nm_col;
        total = 0;
        count = 0;
      }

      void get_sum(int row_in, int col_in,
        const af::const_ref< double, af::c_grid<2> > &data2d,
        const af::const_ref< int, af::c_grid<2> > &mask2d){
          for (int row = row_in - 5; row < row_in + 5; row++){
            for(int col = col_in - 5; col < col_in + 5; col++){
              if(row>=0 and row < max_row and
                 col>=0 and col < max_col ){
                if ( mask2d(row, col) & Background ){
                  total += data2d(row, col);
                  count++;
                }
              }
            }
          }
      }

      double get_avg(){
        double bkgr;
        if ( count > 0) {
          bkgr = total / count;
        } else {
          bkgr = 0;
        }
        return bkgr;
      }


  };


  af::versa< double, af::c_grid<2> > curved_background_flex_2d(
    const af::const_ref< double, af::c_grid<2> > &data2d,
    const af::const_ref< int, af::c_grid<2> > &mask2d) {
        std::size_t ncol=data2d.accessor()[1];
        std::size_t nrow=data2d.accessor()[0];
        af::versa< double, af::c_grid<2> > background2d(data2d.accessor(),0);

        double loc_bkgr = 0;
        int loc_row, loc_col;

        for (int row = 0; row<nrow;row++) {
          for (int col = 0; col<ncol;col++) {
            if ( mask2d(row,col) & Foreground ){

              local_avg lcl_av(nrow, ncol);


              loc_row = nrow - 1;
              loc_col = col;
              lcl_av.get_sum(loc_row, loc_col, data2d, mask2d);


              loc_row = 0;
              loc_col = col;
              lcl_av.get_sum(loc_row, loc_col, data2d, mask2d);


              loc_row = row;
              loc_col = ncol - 1;
              lcl_av.get_sum(loc_row, loc_col, data2d, mask2d);


              loc_row = row;
              loc_col = 0;
              lcl_av.get_sum(loc_row, loc_col, data2d, mask2d);


              loc_bkgr = lcl_av.get_avg();
              background2d(row,col) = loc_bkgr;

            } else {
              background2d(row,col) = data2d(row,col);
            }
          }
        }

    return background2d;
  }

  int get_plane_background_syml_sys_2d(
    const af::const_ref< double, af::c_grid<2> > &data2d,
    const af::const_ref< int, af::c_grid<2> > &mask2d,

     af::ref< double, af::c_grid<2> > mat_a,
     af::ref< double, af::c_grid<2> > vec_b)

  {
        std::size_t ncol=data2d.accessor()[1];
        std::size_t nrow=data2d.accessor()[0];
        af::versa< double, af::c_grid<2> > background2d(data2d.accessor(),0);

        /////////////////////////////////////////////////////////////

        // Counting how many pixels are in Background area
        int counter = 0;
        for (int row = 0; row<nrow;row++) {
          for (int col = 0; col<ncol;col++) {
            if ( mask2d(row, col) & Background ) {
              counter++ ;
            }
          }
        }

        // Building a set of 1D lists with the useful pixels


        double x_lst[counter];
        double y_lst[counter];
        double z_lst[counter];
        counter = 0;
        for (int row = 0; row<nrow; row++) {
          for (int col = 0; col<ncol; col++) {
            if ( mask2d(row, col) & Background ) {

              x_lst[counter] = col + 0.5;
              y_lst[counter] = row + 0.5;
              z_lst[counter] = data2d(row, col);
              counter++ ;
            }
          }
        }
        //elements of the matrix

        double sum_x_sqr = 0.0;
        double sum_x_y   = 0.0;
        double sum_x     = 0.0;
        double sum_y_sqr = 0.0;
        double sum_y     = 0.0;
        double sum_one   = 0.0;

        double sum_x_z   = 0.0;
        double sum_y_z   = 0.0;
        double sum_z     = 0.0;

        //looping trough lists and building the elements of the matrix to solve
        for (int ipos = 0; ipos<counter; ipos++) {
          //std::cout <<"\n"<<ipos<<"  x ="<< x_lst[ipos] << "  y =" << y_lst[ipos] << "  z =" << z_lst[ipos];

          sum_x_sqr += x_lst[ipos] * x_lst[ipos];
          sum_x_y   += x_lst[ipos] * y_lst[ipos];
          sum_x     += x_lst[ipos];
          sum_y_sqr += y_lst[ipos] * y_lst[ipos];
          sum_y     += y_lst[ipos];
          sum_one   += 1.0;
          sum_x_z   += x_lst[ipos] * z_lst[ipos];
          sum_y_z   += y_lst[ipos] * z_lst[ipos];
          sum_z     += z_lst[ipos];
        }
        //std::cout <<"\n";

        mat_a(0,0) = sum_x_sqr;
        mat_a(0,1) = sum_x_y;
        mat_a(0,2) = sum_x;

        mat_a(1,0) = sum_x_y;
        mat_a(1,1) = sum_y_sqr;
        mat_a(1,2) = sum_y;

        mat_a(2,0) = sum_x;
        mat_a(2,1) = sum_y;
        mat_a(2,2) = sum_one;

        vec_b(0,0) = sum_x_z;
        vec_b(0,1) = sum_y_z;
        vec_b(0,2) = sum_z;

        int ok = 0;

    return ok;
  }


  double variance_n_background_from_plane(
    const af::const_ref< double, af::c_grid<2> > &data2d,
    const af::const_ref< int, af::c_grid<2> > &mask2d,
    const af::const_ref< double, af::c_grid<2> > &abc_plane,
    af::ref< double, af::c_grid<2> > background2d){
    //af::versa< double, af::c_grid<2> > background2d(data2d.accessor(),0);
        std::size_t ncol=data2d.accessor()[1];
        std::size_t nrow=data2d.accessor()[0];

        double a = abc_plane(0,0);
        double b = abc_plane(0,1);
        double c = abc_plane(0,2);
        double x, y, z, z_dif, r = 0.0;

        double variance, m = 0, n = 0;

        for (int row = 0; row<nrow;row++) {
          for (int col = 0; col<ncol;col++) {
            x = col + 0.5;
            y = row + 0.5;
            z = a * x + b * y + c;
            if ( mask2d(row,col) & Foreground ){
              background2d(row,col) = z;
              m++;
            } else if ( mask2d(row,col) & Background ){
              background2d(row,col) = data2d(row,col);
              z_dif = z - data2d(row,col);
              r += z_dif * z_dif;
              n++;
            }
          }
        }
        variance = (m / (n-3)) * r;
    return variance;
  }
  /*
   * Illustrative example done by James
   *
   * class InclinedPlaneBackgroundFlex2d {
   * public:
   *   InclinedPlaneBackgroundFlex2d(
   *       const af::const_ref< double, af::c_grid<2> > &data2d,
   *       const af::const_ref< int, af::c_grid<2> > &mask2d,
   *       const af::const_ref< double, af::c_grid<2> > &abc_plane,
   *       af::ref< double, af::c_grid<2> > background2d) {
   *
   *   }
   *
   *   flex_double result() const {
   *     return result_;
   *
   *   }
   *
   *   double sigma() const {
   *     return sigma_;
   *   }
   *
   * };
   */

}}

#endif
