
#ifndef DIALS_REFINEMENT_MAHALANOBIS_H
#define DIALS_REFINEMENT_MAHALANOBIS_H

#include <dials/array_family/scitbx_shared_and_versa.h>
#include <scitbx/array_family/versa_matrix.h>
#include <dials/error.h>

namespace dials { namespace refinement {

  af::shared<double> maha_dist_sq(const af::const_ref<double, af::c_grid<2> > &obs,
                                  const af::const_ref<double> &center,
                                  const af::const_ref<double, af::c_grid<2> > &cov) {
    std::size_t nobs = obs.accessor()[0];
    std::size_t nparam = obs.accessor()[1];

    // check input
    DIALS_ASSERT(nobs > nparam);
    DIALS_ASSERT(center.size() == nparam);
    DIALS_ASSERT(cov.accessor()[0] == nparam);
    DIALS_ASSERT(cov.accessor().is_square());

    // invert covariance matrix
    af::versa<double, af::c_grid<2> > covinv(cov.accessor());
    std::copy(cov.begin(), cov.end(), covinv.begin());
    af::matrix_inversion_in_place(covinv.ref());

    // create output array
    af::shared<double> d2(nobs);

    // loop over observations
    std::size_t k = 0;
    for (std::size_t i = 0; i < nobs; i++) {
      // copy row from obs and shift by means
      af::shared<double> row(nparam);

      for (std::size_t j = 0; j < nparam; j++) {
        row[j] = obs[k] - center[j];
        k++;
      }

      // Mahalanobis distance squared is defined by the matrix product
      //
      // (x - mu)^T [S]^-1 (x - mu)
      //
      // (x - mu) is in 'row' and [S]^-1 is 'covinv'. Perform the right
      // hand part of the product first with matrix_multiply and the
      // second step is just a dot product.
      af::shared<double> prod =
        af::matrix_multiply(covinv.const_ref(), row.const_ref());
      d2[i] = 0.0;
      for (std::size_t j = 0; j < nparam; j++) {
        d2[i] += row[j] * prod[j];
      }
    }
    return d2;
  }

}}  // namespace dials::refinement

#endif  // DIALS_REFINEMENT_MAHALANOBIS_H
