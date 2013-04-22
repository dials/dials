#ifndef DIALS_ALGORITHMS_PEAK_FINDING_LUI_SMOOTHING_H
#define DIALS_ALGORITHMS_PEAK_FINDING_LUI_SMOOTHING_H
#include <iostream>

#include <scitbx/array_family/flex_types.h>

namespace dials { namespace algorithms {
  using scitbx::af::flex_int;

  flex_int smooth_2d(flex_int & data2d, int a) {
        // std::cout << "length =" << data2d.size() << " \n";
        std::size_t ncol=data2d.accessor().all()[1];
        std::size_t nrow=data2d.accessor().all()[0];
        flex_int data2dtmp(data2d);
        flex_int data2dsmoth(data2d.accessor(),0);
        // std::cout <<"a  =" << a << "\n";
        // std::cout <<"ncol =" << ncol << "  nrow =" << nrow <<" \n";
        for (int time = 0; time < a; time++) {
          for (int row = 1; row<nrow-1;row++) {
            for (int col = 1; col<ncol-1;col++) {
                  data2dsmoth(row,col) = (data2dtmp(row - 1, col - 1) + data2dtmp(row - 1, col) + data2dtmp(row - 1, col + 1)
                                + data2dtmp(row  , col - 1) + data2dtmp(row  , col + 1)
                        + data2dtmp(row + 1, col - 1) + data2dtmp(row + 1, col) + data2dtmp(row + 1, col + 1)) / 8.0;
                }
          }
          std::copy(data2dsmoth.begin(),data2dsmoth.end(),data2dtmp.begin() );
        }
    return data2dsmoth;
  }


}}

#endif
//    for times in range(n_times):
//        for row in range(1, n_row - 1, 1):
//            for col in range(1, n_col - 1, 1):
//                data2dsmoth[row, col] = (data2dtmp[row - 1, col - 1] + data2dtmp[row - 1, col] + data2dtmp[row - 1, col + 1]  \
//                                         + data2dtmp[row  , col - 1] + data2dtmp[row  , col + 1]    \
//                                         + data2dtmp[row + 1, col - 1] + data2dtmp[row + 1, col] + data2dtmp[row + 1, col + 1]) / 8.0
//        data2dtmp[:, :] = data2dsmoth[:, :]
