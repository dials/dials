
#ifndef DIALS_ALGORITHMS_STATISTICS_POISSON_TEST_H
#define DIALS_ALGORITHMS_STATISTICS_POISSON_TEST_H

#include <boost/math/distributions.hpp>
#include <dials/error.h>

namespace dials { namespace algorithms {

  inline double poisson_expected_max_counts(double mean, std::size_t nobs) {
    DIALS_ASSERT(nobs > 0);
    boost::math::poisson_distribution<double> d(mean);
    return boost::math::quantile(d, 1 - 1.0 / nobs) + 1;
  }

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_STATISTICS_POISSON_TEST_H
