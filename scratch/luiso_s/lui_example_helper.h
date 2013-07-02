#ifndef DIALS_SCRATCH_LUI_HELPER_H
#define DIALS_SCRATCH_LUI_HELPER_H
#include <iostream>
#include <scitbx/array_family/flex_types.h>
namespace dials { namespace scratch {
  using scitbx::af::flex_int;
  using scitbx::af::flex_double;

  int hello_tst() {
    std::cout << "Hello world (testing)\n";
    int a=5;
    return a;
  }

  flex_int tst_01(flex_int & data2d) {
    std::size_t ncol=data2d.accessor().all()[1];
    std::size_t nrow=data2d.accessor().all()[0];
    flex_int data2dtmp(data2d);
    for (int row = 0; row<=nrow-1;row++) {
      for (int col = 0; col<=ncol-1;col++) {
        data2dtmp(row,col)=row*nrow+col;
      }
    }

    return data2dtmp;

  }
}}

#endif
