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

#include <algorithm>
#include <map>
#include <memory>
#include <vector>
#include <boost/python.hpp>
#include <boost/variant.hpp>
#include <boost/mpl/list.hpp>
#include <boost/mpl/remove_if.hpp>
#include <boost/mpl/transform.hpp>
#include <dials/error.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
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
