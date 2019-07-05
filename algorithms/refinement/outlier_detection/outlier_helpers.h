
#ifndef DIALS_REFINEMENT_OUTLIER_HELPERS_H
#define DIALS_REFINEMENT_OUTLIER_HELPERS_H

#include <dials/error.h>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/gamma.hpp>

namespace dials { namespace refinement {

  double qchisq(double p, double df) {
    // check input is valid
    DIALS_ASSERT(df > 0.0);
    DIALS_ASSERT(p >= 0.0);
    DIALS_ASSERT(p <= 1.0);

    // chi^2 distribution with df d.o.f:
    boost::math::chi_squared_distribution<double> chisq(df);

    // Find x s.t. F(x; df) = p where F is the CDF of the Chi^2 dist with df d.o.f
    return quantile(chisq, p);
  }

  double mcd_consistency(double df, double alpha) {
    // check input is valid
    DIALS_ASSERT(df > 0.0);
    DIALS_ASSERT(alpha >= 0.0);
    DIALS_ASSERT(alpha <= 1.0);

    // qalpha
    double qalpha = qchisq(alpha, df);

    // gamma distribution with shape parameter df/2 + 1
    boost::math::gamma_distribution<> gamma(df / 2 + 1);

    // 1/c_alpha
    double calpha_inv = cdf(gamma, qalpha / 2) / alpha;

    return 1.0 / calpha_inv;
  }
}}  // namespace dials::refinement

#endif  // DIALS_REFINEMENT_OUTLIER_HELPERS_H
