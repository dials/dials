#ifndef DIALS_ALGORITHMS_LUI_2D_BACKGROUND_H
#define DIALS_ALGORITHMS_LUI_2D_BACKGROUND_H
#include <iostream>
#include <scitbx/array_family/flex_types.h>

namespace dials { namespace algorithms {
  using scitbx::af::flex_int;

  flex_int background_subtract_2d(flex_int & data2d, flex_int & mask2d) {
        // std::cout << "length =" << data2d.size() << " \n";
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
            }
          }
        }
    return background2d;
  }


}}
/*


def flat_background_calc_2d(data2d, diffdata2d_ext):

    import numpy
    n_col = numpy.size(data2d[0:1, :])
    n_row = numpy.size(data2d[:, 0:1])
    avg_bkgr_2d = numpy.copy(data2d)

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
            if diffdata2d_ext[row, col] != 0:
                avg_bkgr_2d[row, col] = avg_bkgr


    return avg_bkgr_2d

def curved_background_calc_2d(data2d, mask2d):

    import numpy
    n_col = numpy.size(data2d[0:1, :])
    n_row = numpy.size(data2d[:, 0:1])
    avg_bkgr_2d = numpy.copy(data2d)
    for row in range(0, n_row, 1):
        for col in range(0, n_col, 1):
            if mask2d[row, col] == 1:
                loc_bkgr_tot = 0.0
                loc_bkgr_cont = 0
                if mask2d[n_row - 1, col] == 0:
                    loc_bkgr_tot += data2d[n_row - 1, col]
                    loc_bkgr_cont += 1
                if mask2d[0, col] == 0:
                    loc_bkgr_tot += data2d[0, col]
                    loc_bkgr_cont += 1
                if mask2d[row, n_col - 1] == 0:
                    loc_bkgr_tot += data2d[row, n_col - 1]
                    loc_bkgr_cont += 1
                if mask2d[row, 0] == 0:
                    loc_bkgr_tot += data2d[row, 0]
                    loc_bkgr_cont += 1

                if loc_bkgr_cont > 0:
                    loc_bkgr = loc_bkgr_tot / float(loc_bkgr_cont)
                else:
                    loc_bkgr = 0

                avg_bkgr_2d[row, col] = loc_bkgr


    return avg_bkgr_2d


 */


#endif
