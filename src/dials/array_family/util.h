
#ifndef DIALS_ARRAY_FAMILY_UTIL_H
#define DIALS_ARRAY_FAMILY_UTIL_H

#include <functional>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error>

namespace dials { namespace af {

  /**
   * Copy the indices to std::size_t
   */
  af::shared<std::size_t> int_to_size_t(const const_ref<int> &x) {
    typedef af::shared<std::size_t>::iterator iterator;
    af::shared<std::size_t> result(x.size());
    iterator end =
      std::copy_if(x.begin(), x.end(), result.begin(), std::greater_equal(0));
    DIALS_ASSERT(end == result.end());
    return result;
  }

}}  // namespace dials::af

#endif  // DIALS_ARRAY_FAMILY_UTIL_H
