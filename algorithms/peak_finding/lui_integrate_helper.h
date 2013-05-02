#ifndef DIALS_ALGORITHMS_INTEGRATE_LUI_HELPER_H
#define DIALS_ALGORITHMS_INTEGRATE_LUI_HELPER_H
#include <iostream>
#include <scitbx/array_family/flex_types.h>
namespace dials { namespace algorithms {
  using scitbx::af::flex_int;
  flex_int ref_2d(flex_int & data2d, int tot_times) {
    int ncol=data2d.accessor().all()[1];
    int nrow=data2d.accessor().all()[0];
    flex_int data2dsmoth(data2d.accessor(),0);
    std::cout << "times ="<< tot_times << "\n";
    std::cout << "size(x) =" << ncol <<"\n";
    std::cout << "size(y) =" << nrow <<"\n";

    for (int row = 1; row<nrow-1;row++) {
      for (int col = 1; col<ncol-1;col++) {
        data2dsmoth(row,col) = tot_times;
      }
    }


    return data2dsmoth;
  }

}}

#endif
