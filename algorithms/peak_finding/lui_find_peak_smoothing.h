#ifndef DIALS_ALGORITHMS_PEAK_FINDING_LUI_SMOOTHING_H
#define DIALS_ALGORITHMS_PEAK_FINDING_LUI_SMOOTHING_H
#include <iostream>

#include <scitbx/array_family/flex_types.h>

namespace dials { namespace algorithms {
  using scitbx::af::flex_int;
  using scitbx::af::flex_double;
  flex_double smooth_2d(flex_double & data2d, int tot_times) {
        // std::cout << "length =" << data2d.size() << " \n";
        std::size_t ncol=data2d.accessor().all()[1];
        std::size_t nrow=data2d.accessor().all()[0];
        flex_double data2dtmp(data2d);
        flex_double data2dsmoth(data2d.accessor(),0);
        float tot_i,cont;
        // std::cout <<"times =" << tot_times << "\n";
        // std::cout <<"ncol =" << ncol << "  nrow =" << nrow <<" \n";
        std::cout <<"2D in -> tst" << " \n";

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


  flex_int smooth_3d(flex_int & data3d, int tot_times) {
    std::size_t ncol=data3d.accessor().all()[2];
    std::size_t nrow=data3d.accessor().all()[1];
    std::size_t nfrm=data3d.accessor().all()[0];
    float tot_i,cont;
    // int tot_i,cont;
    flex_int data3dtmp(data3d);
    flex_int data3dsmoth(data3d.accessor(),0);
    // std::cout << "length =" << data3d.size() << " \n";
    // std::cout <<"times =" << tot_times << "\n";
    // std::cout <<"ncol =" << ncol << "\n";
    // std::cout <<"nrow =" << nrow <<" \n";
    // std::cout <<"nfrm =" << nfrm <<" \n";
    std::cout <<"[dimensions] = 3D " << " \n";
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
