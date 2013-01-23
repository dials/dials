
#ifndef DIALS_ARRAY_FAMILY_ARRAY_TYPES_H
#define DIALS_ARRAY_FAMILY_ARRAY_TYPES_H

#include <scitbx/array_family/flex_types.h>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <cctbx/miller.h>

namespace dials { namespace af {

// Copied from scitbx/array_family/flex_types.h to give convenient access to 
// oft used flex array types
#define DIALS_ARRAY_FAMILY_FLEX_TYPEDEFS(T, N) \
    typedef scitbx::af::flex<T >::type flex_ ## N; \
    typedef scitbx::af::flex_const_ref<T >::type flex_ ## N ## _const_ref; \
    typedef scitbx::af::flex_ref<T >::type flex_ ## N ## _ref;

    DIALS_ARRAY_FAMILY_FLEX_TYPEDEFS(scitbx::vec2 <int>, vec2_int)
    DIALS_ARRAY_FAMILY_FLEX_TYPEDEFS(scitbx::vec2 <double>, vec2_double)
    DIALS_ARRAY_FAMILY_FLEX_TYPEDEFS(scitbx::vec3 <int>, vec3_int)
    DIALS_ARRAY_FAMILY_FLEX_TYPEDEFS(scitbx::vec3 <double>, vec3_double)
    DIALS_ARRAY_FAMILY_FLEX_TYPEDEFS(cctbx::miller::index <>, miller_index)

#undef DIALS_ARRAY_FAMILY_FLEX_TYPEDEFS

}} // namespace dials::af

#endif // DIALS_ARRAY_FAMILY_ARRAY_TYPES_H