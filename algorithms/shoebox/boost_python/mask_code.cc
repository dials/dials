/*
 * mask_code.cc
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/algorithms/shoebox/mask_code.h>

namespace dials { namespace algorithms { namespace shoebox { namespace boost_python {

  using namespace boost::python;

  void export_mask_code() {
    enum_<MaskCode>("MaskCode")
      .value("Valid", Valid)
      .value("Background", Background)
      .value("Foreground", Foreground)
      .value("Strong", Strong)
      .value("BackgroundUsed", BackgroundUsed)
      .value("Overlapped", Overlapped);
  }

}}}}  // namespace dials::algorithms::shoebox::boost_python
