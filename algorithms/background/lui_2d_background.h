#ifndef DIALS_ALGORITHMS_LUI_2D_BACKGROUND_H
#define DIALS_ALGORITHMS_LUI_2D_BACKGROUND_H
#include <iostream>
#include <scitbx/array_family/flex_types.h>

namespace dials { namespace algorithms {
  using scitbx::af::flex_int;

  flex_int background_subtract_2d(flex_int & data2d) {
        // std::cout << "length =" << data2d.size() << " \n";
        std::size_t ncol=data2d.accessor().all()[1];
        std::size_t nrow=data2d.accessor().all()[0];
        flex_int data2d_aft(data2d.accessor(),0);

        for (int row = 0; row<nrow;row++) {
          for (int col = 0; col<ncol;col++) {
            data2d_aft(row,col)=data2d(row,col);
          }
        }

        std::cout << "Hi _____________________________ there \n";
    return data2d_aft;
  }


}}
/*
def flat_background_subtraction_2d(data2d, diffdata2d_ext):


    import numpy
    n_col = numpy.size(data2d[0:1, :])
    n_row = numpy.size(data2d[:, 0:1])

    tot_bkgr = 0.0
    cont = 0.0
    for row in range(0, n_row, 1):
        for col in range(0, n_col, 1):
            if diffdata2d_ext[row, col] == 0:
                cont += 1
                tot_bkgr += data2d[row, col]
    if tot_bkgr > 0 and cont > 0:
        avg_bkgr = tot_bkgr / cont
    else:
        avg_bkgr = 0
    #print 'avg_bkgr=', avg_bkgr
    for row in range(0, n_row, 1):
        for col in range(0, n_col, 1):
            if diffdata2d_ext[row, col] == 1 and data2d[row, col] > avg_bkgr:
                data2d[row, col] = data2d[row, col] - avg_bkgr
            else:
                data2d[row, col] = 0

    return avg_bkgr

*/
#endif
