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
        double cont=0.0, tot_bkgr = 0.0;
        double avg_bkgr = 0;
        for (int row = 0; row<nrow;row++) {
          for (int col = 0; col<ncol;col++) {
            if ( mask2d(row,col) & Background ){
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
            if ( mask2d(row,col) & Foreground ){
              background2d(row,col) = avg_bkgr;
            } else {
              background2d(row,col) = data2d(row,col);
            }
          }
        }
    return background2d;
  }
  af::versa< double, af::c_grid<2> > curved_background_flex_2d(
    const af::const_ref< double, af::c_grid<2> > &data2d,
    const af::const_ref< int, af::c_grid<2> > &mask2d) {
        std::size_t ncol=data2d.accessor()[1];
        std::size_t nrow=data2d.accessor()[0];
        af::versa< double, af::c_grid<2> > background2d(data2d.accessor(),0);
        double loc_bkgr_cont, loc_bkgr_tot;
        double loc_bkgr = 0;

        for (int row = 0; row<nrow;row++) {
          for (int col = 0; col<ncol;col++) {
            if ( mask2d(row,col) & Foreground ){
              loc_bkgr_tot = 0.0;
              loc_bkgr_cont = 0.0;
              if ( mask2d(nrow - 1, col) & Background ){
                 loc_bkgr_tot += data2d(nrow - 1, col);
                 loc_bkgr_cont++;
              }
              if ( mask2d(0, col) & Background){
                loc_bkgr_tot += data2d(0, col);
                loc_bkgr_cont++;
              }
              if ( mask2d(row, ncol - 1) & Background ){
                loc_bkgr_tot += data2d(row, ncol - 1);
                loc_bkgr_cont++;
              }
              if ( mask2d(row, 0) & Background) {
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
