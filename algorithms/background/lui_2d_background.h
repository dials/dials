#ifndef DIALS_ALGORITHMS_LUI_2D_BACKGROUND_H
#define DIALS_ALGORITHMS_LUI_2D_BACKGROUND_H
#include <iostream>
#include <scitbx/array_family/flex_types.h>

namespace dials { namespace algorithms {
  using scitbx::af::flex_int;

  flex_int flat_background_flex_2d(flex_int & data2d, flex_int & mask2d) {
        std::size_t ncol=data2d.accessor().all()[1];
        std::size_t nrow=data2d.accessor().all()[0];
        flex_int background2d(data2d.accessor(),0);
        float cont=0.0, tot_bkgr = 0.0;
        int avg_bkgr = 0;
        for (int row = 0; row<nrow;row++) {
          for (int col = 0; col<ncol;col++) {
            if ( mask2d(row,col) == 0 ){
              tot_bkgr += data2d(row,col);
              cont++;
            }
            if ( tot_bkgr > 0 and cont > 0) {
                avg_bkgr = tot_bkgr / cont;
            }
          }
        }
        for (int row = 0; row<nrow;row++) {
          for (int col = 0; col<ncol;col++) {
            if ( mask2d(row,col) == 1 ){
              background2d(row,col) = avg_bkgr;
            } else {
              background2d(row,col) = data2d(row,col);
            }
          }
        }
    return background2d;
  }
  flex_int curved_background_flex_2d(flex_int & data2d, flex_int & mask2d) {
        std::size_t ncol=data2d.accessor().all()[1];
        std::size_t nrow=data2d.accessor().all()[0];
        flex_int background2d(data2d.accessor(),0);
        float loc_bkgr_cont, loc_bkgr_tot;
        int loc_bkgr = 0;

        for (int row = 0; row<nrow;row++) {
          for (int col = 0; col<ncol;col++) {
            if ( mask2d(row,col) == 1 ){
              loc_bkgr_tot = 0.0;
              loc_bkgr_cont = 0.0;
              if ( mask2d(nrow - 1, col) == 0 ){
                 loc_bkgr_tot += data2d(nrow - 1, col);
                 loc_bkgr_cont++;
              }
              if ( mask2d(0, col) == 0){
                loc_bkgr_tot += data2d(0, col);
                loc_bkgr_cont++;
              }
              if ( mask2d(row, ncol - 1) == 0 ){
                loc_bkgr_tot += data2d(row, ncol - 1);
                loc_bkgr_cont++;
              }
              if ( mask2d(row, 0) == 0) {
                loc_bkgr_tot += data2d(row, 0);
                loc_bkgr_cont++;
              }
              if ( loc_bkgr_cont > 0) {
                loc_bkgr = loc_bkgr_tot / loc_bkgr_cont;
              } else {
                loc_bkgr = 0;
              }
              background2d(row,col) = loc_bkgr;
            } else {
              background2d(row,col) = data2d(row,col);
            }
          }
        }
    return background2d;
  }

}}

#endif
