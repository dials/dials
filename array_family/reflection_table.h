/*
 * reflection_table.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ARRAY_FAMILY_REFLECTION_TABLE_H
#define DIALS_ARRAY_FAMILY_REFLECTION_TABLE_H

#include <dials/array_family/flex_table.h>
#include <dials/model/data/shoebox.h>
#include <dials/model/data/transformed_shoebox.h>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/mat3.h>
#include <scitbx/vec3.h>
#include <scitbx/vec2.h>
#include <cctbx/miller.h>

namespace dials { namespace af {

  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::mat3;
  using scitbx::af::int6;
  using model::Shoebox;
  using model::TransformedShoebox;

  typedef flex_type_generator<
    bool,
    int,
    std::size_t,
    double,
    std::string,
    vec2<double>,
    vec3<double>,
    mat3<double>,
    int6,
    cctbx::miller::index<>,
    Shoebox<>,
    TransformedShoebox
  >::type reflection_table_types;

  typedef flex_table<reflection_table_types> reflection_table;

  enum Flags {
    Predicted        = (1 << 0),
    Observed         = (1 << 1),
    Indexed          = (1 << 2),
    UsedInRefinement = (1 << 3),
    ReferenceSpot    = (1 << 4),
    Integrated       = (1 << 5),
    DontProcess      = (1 << 6),
  };

}} // namespace dials::af

#endif // DIALS_ARRAY_FAMILY_REFLECTION_TABLE_H
