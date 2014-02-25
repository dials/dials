#ifndef DIALS_ALGORITHMS_PEAK_FINDING_LUI_HELPER_H
#define DIALS_ALGORITHMS_PEAK_FINDING_LUI_HELPER_H
#include <iostream>
#include <dials/array_family/scitbx_shared_and_versa.h>
namespace dials { namespace algorithms {

  af::versa< double, af::c_grid<2> > find_mask_2d(
        const af::const_ref< double, af::c_grid<2> > &data2d
      , const af::const_ref< double, af::c_grid<2> > &data2dsmoth
      , int tot_times) {
        af::versa< double, af::c_grid<2> > mask2d(data2d.accessor(),0);
        std::size_t ncol=data2d.accessor()[1];
        std::size_t nrow=data2d.accessor()[0];
        std::cout << "finding mask 01 \n";
        for (int row = 1; row<nrow-1;row++) {
          for (int col = 1; col<ncol-1;col++){
            if (data2d(row,col) > data2dsmoth(row,col) + 5.0){
              mask2d(row,col) = 1;
            }else{
              mask2d(row,col) = 0;
            }
          }
        }

    return mask2d;
  }

}}

#endif
