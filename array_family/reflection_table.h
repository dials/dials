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
    Shoebox<>
  >::type reflection_table_types;

  typedef flex_table<reflection_table_types> reflection_table;

  enum Flags {
    Predicted        = (1 << 0),
    Observed         = (1 << 1),
    Indexed          = (1 << 2),
    UsedInRefinement = (1 << 3),
    Strong           = (1 << 4),

    // Role in integration
    ReferenceSpot    = (1 << 5),
    DontIntegrate    = (1 << 6),

    // Integated
    IntegratedSum    = (1 << 7),
    IntegratedPrf    = (1 << 8),
    Integrated       = IntegratedSum | IntegratedPrf,

    // Bad shoebox
    Overloaded       = (1 << 9),
    OverlappedBg     = (1 << 10),
    OverlappedFg     = (1 << 11),
    InPowderRing     = (1 << 12),
    BadShoebox       = Overloaded | OverlappedBg | OverlappedFg | InPowderRing,

    // Bad background
    LargeBgVariation = (1 << 13),
    LargeBgGradient  = (1 << 14),
    LargeBgLevel     = (1 << 15),
    BadBackground    = LargeBgVariation | LargeBgGradient | LargeBgLevel,

    // Bad shape
    LargePosDiff     = (1 << 16),
    BadShape         = LargePosDiff,

    // Bad integrated intensity
    LargeSumPrfDiff  = (1 << 17),
    NegativeSum      = (1 << 18),
    NegativePrf      = (1 << 19),
    PoorProfileFit   = (1 << 20),
    BadIntegration   = LargeSumPrfDiff | NegativeSum | NegativePrf | PoorProfileFit,

    // Bad spot
    BadSpot = BadShoebox | BadBackground | BadShape | BadIntegration
  };

}} // namespace dials::af

#endif // DIALS_ARRAY_FAMILY_REFLECTION_TABLE_H
