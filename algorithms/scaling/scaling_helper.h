#ifndef DIALS_SCRATCH_SCALING_SCALING_HELPER_H
#define DIALS_SCRATCH_SCALING_SCALING_HELPER_H

#include <scitbx/sparse/matrix.h>

namespace dials_scratch { namespace scaling {
  
  /**
   * Elementwise squaring of a matrix
   */
  scitbx::sparse::matrix<double> elementwise_square(scitbx::sparse::matrix<double> m){

    scitbx::sparse::matrix<double> result(m.n_rows(), m.n_cols());

    // outer loop iterate over the columns
    for (std::size_t j=0; j < m.n_cols(); j++) {

      // inner loop iterate over the non-zero elements of the column
      for (scitbx::sparse::matrix<double>::row_iterator p=m.col(j).begin(); p != m.col(j).end(); ++p)
      {
        std::size_t i = p.index();
        result(i, j) = *p * *p;
      }
    }
    return result;
  }

}} // dials_scratch::scaling

#endif // DIALS_SCRATCH_SCALING_SCALING_HELPER_H
