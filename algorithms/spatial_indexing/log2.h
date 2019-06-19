#ifndef DIALS_ALGORITHMS_SPATIAL_INDEXING_LOG2
#define DIALS_ALGORITHMS_SPATIAL_INDEXING_LOG2

#include <cmath>

namespace dials { namespace algorithms {

// FIXME Windows doesn't have a log2. Couldn't find one anywhere else but should
// probably move this to a better location.
#ifdef _WIN32
  inline double log2(double x) {
    return log(x) / log(2.0);
  }
#endif

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_SPATIAL_INDEXING_LOG2
