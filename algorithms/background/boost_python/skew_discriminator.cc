/*
 * skew_discriminator.cc
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
#include <dials/algorithms/background/skew_discriminator.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_skew_discriminator()
  {
    class_<SkewDiscriminator, bases<DiscriminatorStrategy> >(
        "SkewDiscriminator")
      .def(init<>())
      .def("__call__", &SkewDiscriminator::operator());
  }

}}} // namespace = dials::algorithms::boost_python
