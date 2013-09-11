#ifndef DIALS_ALGORITHMS_PEAK_FINDING_LUI_SMOOTHING_H
#define DIALS_ALGORITHMS_PEAK_FINDING_LUI_SMOOTHING_H
#include <iostream>

#include <dials/array_family/scitbx_shared_and_versa.h>

namespace dials { namespace algorithms {

  af::versa< double, af::c_grid<2> > smooth_2d(
      const af::const_ref< double, af::c_grid<2> > &data2d, int tot_times) {

    std::size_t ncol=data2d.accessor()[1];
        std::size_t nrow=data2d.accessor()[0];
        af::versa< double, af::c_grid<2> > data2dtmp(data2d.accessor());
        af::versa< double, af::c_grid<2> > data2dsmoth(data2d.accessor(),0);
        long double tot_i,cont;

        for (int i = 0; i < data2d.size(); ++i) {
          data2dtmp[i] = data2d[i];
        }

        // scanning trough all pixels
        for (int time = 0; time < tot_times; time++) {
          for (int row = 1; row<nrow-1;row++) {
            for (int col = 1; col<ncol-1;col++) {

              // average of all surrounding pixels
              tot_i=0.0;
              cont=0.0;
              for (int lp_row = row - 1; lp_row <= row + 1; lp_row++) {
                for (int lp_col = col - 1; lp_col <= col + 1; lp_col++) {
                  if (lp_row != row or lp_col != col){
                    tot_i += data2dtmp(lp_row, lp_col);
                    cont++;
                  }
                }
              }
              data2dsmoth(row,col) = tot_i/cont;
            }
          }
          std::copy(data2dsmoth.begin(),data2dsmoth.end(),data2dtmp.begin() );
        }

    return data2dsmoth;
  }


  af::versa< double, af::c_grid<3> > smooth_3d(
      const af::const_ref< double, af::c_grid<3> > & data3d, int tot_times) {
    std::size_t ncol=data3d.accessor()[2];
    std::size_t nrow=data3d.accessor()[1];
    std::size_t nfrm=data3d.accessor()[0];
    float tot_i,cont;
    af::versa< double, af::c_grid<3> > data3dtmp(data3d.accessor());
    af::versa< double, af::c_grid<3> > data3dsmoth(data3d.accessor(),0);
    for (int i = 0; i < data3d.size(); ++i) {
      data3dtmp[i] = data3d[i];
    }

                                          // scanning the block
    for (int time = 0; time < tot_times; time++) {
      for (int frm = 1; frm<nfrm-1;frm++) {
        for (int row = 1; row<nrow-1;row++) {
          for (int col = 1; col<ncol-1;col++) {

                                          // average of all surrounding pixels
            tot_i=0.0;
            cont=0;
            for (int lp_frm = frm - 1; lp_frm <= frm + 1; lp_frm++) {
              for (int lp_row = row - 1; lp_row <= row + 1; lp_row++) {
                for (int lp_col = col - 1; lp_col <= col + 1; lp_col++) {
                  if (lp_frm != frm or lp_row != row or lp_col != col){
                    tot_i += data3dtmp(lp_frm, lp_row, lp_col);
                    cont++;
                  }
                }
              }
            }
            data3dsmoth(frm,row,col) = tot_i/cont;
          }
        }
      }
      std::copy(data3dsmoth.begin(),data3dsmoth.end(),data3dtmp.begin() );
    }
    return data3dsmoth;
  }


}}

#endif
