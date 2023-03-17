/*
 * flex_table.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ARRAY_FAMILY_FLEX_TABLE_H
#define DIALS_ARRAY_FAMILY_FLEX_TABLE_H

#include <dxtbx/array_family/flex_table.h>

namespace dials { namespace af {

  using flex_table
    [[deprecated("flex_table has moved to dxtbx/array_family/flex_table.h")]] =
      dxtbx::af::flex_table;
  using flex_type_generator
    [[deprecated("flex_type_generator has moved to dxtbx/array_family/flex_table.h")]] =
      dxtbx::af::flex_type_generator;

}}  // namespace dials::af

#endif  // DIALS_ARRAY_FAMILY_FLEX_TABLE_H
