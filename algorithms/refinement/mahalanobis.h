
#ifndef DIALS_REFINEMENT_MAHALANOBIS_H
#define DIALS_REFINEMENT_MAHALANOBIS_H

#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace refinement {

  af::shared<double> maha_dist_sq(
      const af::const_ref< double, af::c_grid<2> > &obs,
      const af::const_ref< double> &center,
      const af::const_ref< double, af::c_grid<2> > &cov){

    int nobs = obs.accessor()[0];
    int nparam = obs.accessor()[1];

    //std::cout << "nobs " << nobs << "\n";
    //std::cout << "nparam " << nparam << "\n";

    // check input
    DIALS_ASSERT (nobs > nparam);
    DIALS_ASSERT (center.size() == nparam);
    DIALS_ASSERT (cov.accessor()[0] == nparam);
    DIALS_ASSERT (cov.accessor()[1] == nparam);

    af::shared<double> d2(nobs);

    for(std::size_t i = 0; i < nobs; i++) {
       d2[i] = 0.0;
    }
    return d2;
  }

}} // namespace dials::refinement

#endif // DIALS_REFINEMENT_MAHALANOBIS_H
