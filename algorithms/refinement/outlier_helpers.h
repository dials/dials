
#ifndef DIALS_REFINEMENT_OUTLIER_HELPERS_H
#define DIALS_REFINEMENT_OUTLIER_HELPERS_H

#include <dials/error.h>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/gamma.hpp>

namespace dials { namespace refinement {

  double mcd_consistency(double p,
                         double alpha){

    // check input is valid
    DIALS_ASSERT (p > 0.0);
    DIALS_ASSERT (alpha >= 0.0);
    DIALS_ASSERT (alpha <= 1.0);

    // chi^2 distribution with p d.o.f:
    boost::math::chi_squared_distribution<double> chisq(p);

    // qalpha
    double qalpha = quantile(chisq, alpha);

    // gamma distribution with shape parameter p/2 + 1
    boost::math::gamma_distribution<> gamma(p/2 + 1);

    // 1/c_alpha
    double calpha_inv = cdf(gamma, qalpha/2) / alpha;

    return 1.0/calpha_inv;
  }
}} // namespace dials::refinement

#endif // DIALS_REFINEMENT_OUTLIER_HELPERS_H
