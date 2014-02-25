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
        std::cout << "finding mask 01 \n";
    return mask2d;
  }

}}

#endif
