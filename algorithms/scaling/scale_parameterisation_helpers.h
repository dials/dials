/*
 * scale_parameterisation_helpers.h
 *
 *  Copyright (C) (2016) STFC Rutherford Appleton Laboratory, UK.
 *
 *  Author: David Waterman
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef DIALS_ALGORITHMS_SCALING_SCALE_PARAMETERISATION_HELPERS_H
#define DIALS_ALGORITHMS_SCALING_SCALE_PARAMETERISATION_HELPERS_H

#include <dials/array_family/scitbx_shared_and_versa.h>
#include <scitbx/sparse/matrix.h>
#include <dials/error.h>

namespace dials { namespace scaling {

  /**
   * Multiply the rows of a sparse matrix by the elements of a vector.
   * @param m The double precision sparse matrix
   * @param v The vector of values to scale each row by
   */
  scitbx::sparse::matrix<double> row_multiply(scitbx::sparse::matrix<double> m,
                                              af::const_ref<double> v){

    DIALS_ASSERT(m.n_rows() == v.size());

    // call compact to ensure that each elt of the matrix is only defined once
    m.compact();

    scitbx::sparse::matrix<double> result(m.n_rows(), m.n_cols());

    // outer loop iterate over the columns
    for (std::size_t j=0; j < m.n_cols(); j++) {

      // inner loop iterate over the non-zero elements of the column
      for (scitbx::sparse::matrix<double>::row_iterator p=m.col(j).begin(); p != m.col(j).end(); ++p)
      {
        std::size_t i = p.index();
        result(i, j) = *p * v[i];
      }
    }
    return result;
  }

}} // namespace dials::scaling

#endif // DIALS_ALGORITHMS_SCALING_SCALE_PARAMETERISATION_HELPERS_H
