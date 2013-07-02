#ifndef DIALS_SCRATCH_LUI_HELPER_H
#define DIALS_SCRATCH_LUI_HELPER_H
#include <iostream>
#include <scitbx/vec2.h>
#include <scitbx/array_family/flex_types.h>
namespace dials { namespace scratch {
  using scitbx::af::flex_int;
  using scitbx::af::flex_double;
  using scitbx::vec2;

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

  vec2<double> raw_2d_cut(flex_double & data2d, flex_int & mask2d,
      flex_double & background2d) {

      //std::size_t ncol=data2d.accessor().all()[1];
      //std::size_t nrow=data2d.accessor().all()[0];
      vec2<double> integr_data(0,1);
      /*
      for (int row = 0; row<=nrow-1;row++) {
        for (int col = 0; col<=ncol-1;col++) {
          data2dtmp(row,col)=row*nrow+col;

        }
      }
      */
      integr_data[0]=5;   // intensity
      integr_data[1]=15;  // intensity variance
      return integr_data;

    }


}}
/*

    import numpy
    n_col = numpy.size(data2d[0:1, :])
    n_row = numpy.size(data2d[:, 0:1])

    #print data2d
    i_tot = 0
    npix_bkgr = 0
    npix_mask = 0
    cont = 0
    tot_bkgr = 0
    for col in range(n_col):
        for row in range(n_row):

            if mask2d[row, col] == 1 :
                i_tot = i_tot + data2d[row, col] - bkgr2d[row, col]
                npix_mask += 1
            else:
                npix_bkgr += 1

            cont += 1
            tot_bkgr += bkgr2d[row, col]
    if tot_bkgr > 0 and cont > 0 and npix_mask > 0 and npix_bkgr > 0 :  #fix me
        bkgr = tot_bkgr / cont
        sig = numpy.sqrt(i_tot + (1.0 + (npix_mask) / (npix_bkgr)) * (npix_mask * bkgr))
    else:
        bkgr = 0
        sig = numpy.sqrt(i_tot)

    return i_tot, sig

*/

#endif
